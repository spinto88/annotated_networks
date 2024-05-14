CFLAGS = -O
CC = gcc
LIBS = -lgsl -lgslcblas -lm

metadata: readgml.o
	$(CC) $(CFLAGS) -o metadata.e metadata.c $(LIBS) readgml.o

readgml.o: readgml.c network.h Makefile
	$(CC) $(CFLAGS) -c readgml.c
