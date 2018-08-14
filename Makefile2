CFLAGS = -O -g
CC = gcc
LIBS = -lm -lgsl -lgslcblas
TARGET = main

SRC = $(wildcard *.c)
OBJ = $(patsubst %.c, %.o, $(SRC))

.PHONY: all object clean

default: all

object: $(OBJECTS)

serial: $(TARGET)

all: object serial

$(TARGET): $(OBJ)
	$(CC) $^ -o $(TARGET) $(LDFLAGS) $(LIBS)

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@ 

clean:
	rm -rfv *.o
