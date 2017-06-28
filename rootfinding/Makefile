CFLAGS = -Wall -std=gnu99 -O `gsl-config --cflags`
LDLIBS = $$(gsl-config --libs)
CC = clang

all:	data
	cat data

data: main Makefile
		./main > data

main: main.o root.o

.PHONEY: clean
clean: 
		$(RM) main *.o *.txt data *.png
