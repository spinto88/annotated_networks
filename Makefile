CFLAGS = -O -g -fPIC
CC = gcc
LIBS = -lm
TARGET = libc.so

SRC = $(wildcard *.c)
OBJ = $(patsubst %.c, %.o, $(SRC))

.PHONY: all object clean

default: all

object: $(OBJECTS)

serial: $(TARGET)

all: object serial 

$(TARGET): $(OBJ)
	$(CC) -shared $^ -o $(TARGET) $(LDFLAGS) $(LIBS)

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@ 

clean:
	rm -rfv *.o *.so
