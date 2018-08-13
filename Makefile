CFLAGS = -O
CC = gcc
LIBS = -lm -lgsl -lgslcblas

metadata:metadata.o readgml.o
	$(CC) $(CFLAGS) metadata.o readgml.o -o metadata $(LIBS)

metadata.o:metadata.c readgml.h network.h Makefile
	$(CC) $(CFLAGS) -c metadata.c 

readgml.o: readgml.c network.h
	$(CC) $(CFLAGS) -c readgml.c

clean:
	rm -rfv readgml.o metadata.o metadata
