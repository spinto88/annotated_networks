CFLAGS = -O
CC = gcc
LIBS = -lm -lgsl

metadata:metadata.o readgml.o
	$(CC) $(CFLAGS) -o metadata $(LIBS) metadata.o readgml.o

metadata.o:metadata.c readgml.h network.h Makefile
	$(CC) $(CFLAGS) -c metadata.c

readgml.o:readgml.c network.h Makefile
	$(CC) $(CFLAGS) -c readgml.c
