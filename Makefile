CC = gcc

# -static not working on Mac
CFLAGS := $(CFLAGS) -g -pipe -Wall -ansi -O2 
# := expands variables (support for environment variables)
INCLUDES := -I./vienna-rna-include
LDFLAGS := $(LDFLAGS) -lRNA -lsquid -lm 


TARGET = scif
SOURCES=scif.c alifold.c
OBJECTS=scif.o alifold.o

$(TARGET):
	$(CC) $(CFLAGS) $(INCLUDES) $(SOURCES) $(LDFLAGS) -o $(TARGET)

clean:
	rm -f $(TARGET) $(OBJECTS)
