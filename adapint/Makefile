CFLAGS = -Wall -std=gnu99 -O `gsl-config --cflags`
LDLIBS = $$(gsl-config --libs)
CC = gcc

all:  data
	cat data

data: main
	./main 1> data


main: main.o adaptint.o

.PHONEY: clean
clean: 
	$(RM) main *.o *.txt data *.png
