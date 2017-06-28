CFLAGS = -Wall -std=gnu99 -O `gsl-config --cflags`
LDLIBS = $$(gsl-config --libs)
CC = gcc

all: data plot.png

plot.png: plot.gpi data
	gnuplot plot.gpi

data: main
	./main 1> data


main: main.o montecarlo.o

.PHONEY: clean
clean: 
	$(RM) main *.o *.txt data *.png
