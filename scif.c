/***
 *
 * -- scif
 *
 * Compute the structure conservation index fast
 *
 *
 * Copyright (C) 2005 Andreas Wilm <wilm@biophys.uni-duesseldorf.de>
 *
 * This file is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 ***/

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
/*#include <ctype.h>*/
#include <stdarg.h>
#include <getopt.h>

#include <squid.h>
#include "AliFold/alifold.h"
#include "fold.h"



#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif


/* TMP
  viennarnadir=/home/opal/wilm/debian/vienna-rna/vienna-rna-1.5/;
  gcc -o scif -Wall -ansi scif.c alifold.c \
     -I$viennarnadir -I$viennarnadir/H \
     -L/usr/lib \
     -lsquid -lRNA -lm
*/



/* print statements: taken from squicl-0.2.8 (squicl/src/utils.h)
 *  append trailing "\n"
 *  use at least one fmt+string
 */
int debug=0;
int be_verbose=0;

/* print only if debug is true*/
#define DEBUG_P(fmt, args...)     printk(stdout, debug, "DEBUG(%s|%s): " fmt, __FILE__, __FUNCTION__, ## args)
/* print only if verbose is true*/
#define VERBOSE_P(fmt, args...)   printk(stdout, be_verbose, fmt, ## args)
/* always warn to stderr */
#define WARN_P(fmt, args...)      printk(stdout, 1, "WARNING(%s): " fmt, __FILE__, ## args)
/* always print errors to stderr*/
#define ERROR_P(fmt, args...)     printk(stderr, 1, "ERROR(%s|%s): " fmt, __FILE__, __FUNCTION__, ## args)
/* always print temp debugging */
#define TMPDEBUG_P(fmt, args...)  printk(stdout, 1, "TMPDEBUG(%s|%s|%d): " fmt, __FILE__, __FUNCTION__, __LINE__, ## args)
/* always print fixme's */
#define FIXME_P(fmt, args...)  printk(stdout, 1, "FIXME(%s|%s|%d): " fmt, __FILE__, __FUNCTION__, __LINE__, ## args)


#define MYNAME "scif"


void
Usage(char *errmsg);
int
printk(FILE *stream, int flag, const char *fmt, ...);
void
MsaToCase(MSA *msa, char c);
float
MSAlifold(MSA *msa, char *structure);
float
Mfe(char *seq, char *structure);
char *
Dealign(const char *aln_seq);
MSA *
MsaRead( char *afile);




/* external functions:
 *
 * squid:isgap
 * squid:MSAFileOpen
 * squid:MSAFileRead
 * squid:MSAFileClose
 * squid:ToRNA
 * squid:MSAMingap
 * squid:FileExists
 * squid:s2upper
 * "vienna-rna:alifold"
 * vienna-rna:fold
 * vienna-rnafree_arrays
 *
 */





/***   Usage   ***
 *
 * FIXME:description
 *
 * IN:
 * OUT:
 * SIDEEFFECTS:
 * NOTES:
 */
void
Usage(char *errmsg)
{
    if (errmsg!=NULL)
        fprintf(stderr, "ERROR: %s\n", errmsg);
    fprintf(stderr, "%s: fast computation of the tructure Conservation Index\n\n", MYNAME);
    fprintf(stderr, "  see Washietl 2005 PNAS: Fast and reliable prediction of noncoding RNAs:\n");
    fprintf(stderr, "  rnaalifold consensus mfe / avg(single rnafold mfe)\n");
    fprintf(stderr, "\nUsage:\n [option] %s [msa-file]\n", MYNAME);
    fprintf(stderr, "Options: -d         : debug\n");
    fprintf(stderr, "         -v         : be verbose\n");
    fprintf(stderr, "         -c <float> : set covariance term for alifold (defaults to 1)\n");
}
/***   Usage   ***/



/***   printk
 *
 * Taken from the Linux kernel source
 * and slightly modified
 *
 * flag: boolean: print or don't
 *       
 */
int
printk(FILE *stream, int bool_flag, const char *fmt, ...)
{                
    va_list args;
    static char printk_buf[8192];
    int printed_len=0;
 
    if (bool_flag) {
        /* Emit the output into the temporary buffer */
        va_start(args, fmt);
        printed_len = vsnprintf(printk_buf, sizeof(printk_buf), fmt, args);
        va_end(args);

        fprintf(stream, "%s", printk_buf);
        fflush(stream);        
    }
    return printed_len;
}
/***   printk   ***/



/***   MsaToCase   ***
 *
 * IN:
 * OUT:
 * SIDEEFFECTS:
 * NOTES:
 *  Taken from squicl-0.2.8  (squicl/src/msa.c)
 */
void
MsaToCase(MSA *msa, char c)
{
    int i;
    
    if (c!='l' && c!='u') {
        fprintf(stderr, "Invalid case %c as argument\n.", c);
        return;
    }
    if (c!='l' && c!='u') {
        WARN_P("Ignoring invalid case %c as argument\n.", c);
        return;
    }

    for (i=0; i<msa->nseq; i++) {
        if (c=='l')
            s2lower(msa->aseq[i]);
        else if (c=='u')
            s2upper(msa->aseq[i]);
    }
}
/***   MsaToCase   ***/



/***   MSAlifold   ***
 *
 * Just a wrapper to alifold
 *
 * IN:
 *  MSA alignment to fold
 *  struct: preallocated string for dotbracket structure notation
 * OUT:
 *   energy in kcal/mol
 * SIDEEFFECTS:
 * NOTES:
 */
float
MSAlifold(MSA *msa, char *structure)
{
    float energy;
    char **aln;
    int i;
    
    if (msa==NULL || structure==NULL)
        return 0.0;

    /* init is called by fold itself, in case seqleng>initlen
       init_alifold(seqlen);
    */
    
    /* BEWARE: alifold determines nseq by trailing NULL pointer ! */
    aln = (char**) MallocOrDie((msa->nseq+1)*sizeof(char*));
    for (i = 0; i < msa->nseq; i++) {
        aln[i] = (char*) MallocOrDie(1*sizeof(char**));
        aln[i] = msa->aseq[i];
    }
    aln[i] = NULL; /* see note above */
    
    energy = alifold(aln, structure);
    DEBUG_P("energy = %f\n", energy);

    /* don't call: free_alifold_arrays() maybe we come back later */
    free(aln);

    return energy;
}
/***   MSAlifold   ***/



/***   Mfe   ***
 *
 * Just a wrapper to RNAfold's fold
 *
 * IN:
 *  seq: (dealigned) rna sequence to fold
 *  struct: preallocated string for dotbracket structure notation
 * OUT:
 *   energy in kcal/mol
 * SIDEEFFECTS:
 * NOTES:
 *   Taken from squicl-0.2.8  (squicl/src/rnafold.c)
 */
float
Mfe(char *seq, char *structure)
{
    int seqlen;
    float energy;

    if (seq==NULL || structure==NULL)
        return 0.0;
    
    seqlen = strlen(seq);

    /* init is called by fold itself, in case seqleng>initlen
       initialize_fold(seqlen);
    */
    
    /* folds the sequence and returns the minimum free energy in
       kcal/mol the mfe structure in bracket notation (see section
       Representations of Secondary Structures) is returned in
       structure. Sufficient space for string of the same length as
       sequence must be allocated for structure before calling
       fold(). */
    energy = fold(seq, structure);
    DEBUG_P("energy = %f\n", energy);
    free_arrays();
        
    /* don't call: free_arrays() maybe we come back later */
    return energy;
}
/***   Mfe   ***/





/***   Dealign   ***
 *
 * Dealigns provided sequence string and returns the string
 *
 * IN:
 * OUT:
 * SIDEEFFECTS:
 * NOTES:
 *  Caller must free
 *  Taken from squicl-0.2.8  (squicl/src/utils.c)
 */
char *
Dealign(const char *aln_seq)
{
    int  i, counter;
    char *ret, *tmp;

    tmp = (char*) MallocOrDie((strlen(aln_seq)+1) * sizeof(char));
    
    counter=0;
    for (i=0; i<strlen(aln_seq); i++) {
        if (! isgap(aln_seq[i])) {
            tmp[counter] = aln_seq[i];
            counter++;
        }
    }
    tmp[counter] = '\0';

    ret = (char*) strdup(tmp);
    free(tmp);
    DEBUG_P("ret = %s\n", ret);

    return ret;
}
/***   Dealign   ***/





/***   MsaRead ***
 *
 * IN:
 *  afile : name of aligned rna sequence file
 *
 + OUT;
 *  Returns first read MSA* or NULL on error
 *
 * SIDEEFFECTS:
 * NOTES:
 *  Caller must free (MSAFree(msa);)
 *  Taken from squicl-0.2.8  (squicl/src/msa.c)
 *
 */
MSA *
MsaRead(char *afile)
{
    MSAFILE  *afp; /* pointer to open alignment file */
    MSA      *msa; /* multiple sequence alignment    */  
    int       fmt; /* format of afile */

    fmt = SQFILE_UNKNOWN;
    if ((afp = MSAFileOpen(afile, fmt, NULL)) == NULL) {
        fprintf(stderr, "Alignment file %s could not be opened for reading", afile);
        return NULL;
    }
    /* read only first MSA */
    msa = MSAFileRead(afp);
    MSAFileClose(afp);
    
    return msa;
}
/***   MsaRead   ***/





/***   main   ***
 *
 * FIXME:description
 *
 * IN:
 * OUT:
 * SIDEEFFECTS:
 * NOTES:
 */
int
main(int argc, char **argv)
{
    MSA *msa;
    char *fmsa, *dummystruct;
    float mfesum, consmfe;
    float sci;
    int i, c;
    

    while((c=getopt(argc, argv, "vdc:"))!= EOF) {
        switch(c) {
        case 'v':
            be_verbose=1;
            break;
        case 'd':
            debug=1;
            break;
        case 'c':
            i=sscanf(optarg, "%lf", &cv_fact);
            if (!i) {
                Usage("cannot parse covariance factor");
                exit(EXIT_FAILURE);
            }
            break;
        default:
            /*Usage();
              exit(EXIT_FAILURE);*/
            break;
        }
    }

    if (optind+1 == argc) {
        fmsa=argv[optind];
    } else {
        Usage("missing file or non-option argument detected");
        exit(EXIT_FAILURE);
    }


    if ( ! FileExists(fmsa)) {
        char msg[1024];
        sprintf(msg, "Non existant file %s\n", fmsa);
        Usage(msg);
        exit(EXIT_FAILURE);
    }
    msa = MsaRead(fmsa);

    dummystruct = (char *) MallocOrDie((msa->alen+1)*sizeof(char));

    
    /* clean the alignments
     *
     */
    for (i = 0; i < msa->nseq; i++)
        ToRNA(msa->aseq[i]);
    MSAMingap(msa); /* del gap only columns */
    /* squicl::seq::ToStrippedDownIupac $msa ? */
    /* to upper case: vienna-rna folding routines prerequsite! */
    MsaToCase(msa, 'u');

        
    /* get energies for single sequences
     */
    mfesum=0.0;
    for (i=0; i<msa->nseq; i++) {
        char *nalseq;
        float mfe;
        
        nalseq = Dealign(msa->aseq[i]);
        mfe = Mfe(nalseq, dummystruct);
        DEBUG_P("mfe for seq %d %s = %f = %s\n", i, msa->sqname[i], mfe, nalseq);
        mfesum += mfe;

        free(nalseq);
    }
    DEBUG_P("mfesum = %f\n", mfesum);
    
    /* get consensus energy (alifold)
     */
    consmfe = MSAlifold(msa, dummystruct);
    DEBUG_P("consmfe = %f\n", consmfe);

    
    if ((mfesum/msa->nseq)!=0.0) 
        sci = consmfe / (mfesum/msa->nseq);
    else
        sci=0.0;
    fprintf(stdout, "%0.2f\n", sci);

    
    MSAFree(msa);
    free(dummystruct);

    return EXIT_SUCCESS;
}
/***   main   ***/
