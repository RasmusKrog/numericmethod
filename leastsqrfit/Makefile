CFLAGS = -Wall -std=gnu99 -O `gsl-config --cflags`
LDLIBS = $$(gsl-config --libs)
CC = clang
PLOTTER = gnuplot

all: plot.png

plot.png: plot.gpi data
	gnuplot plot.gpi

data: main Makefile
	./main > data

main: main.o linearequation.o leastsquare.o

leastsquares: leastsquare.o linearequation.o

.PHONEY: clean
clean: 
	$(RM) main *.o *.txt data *.png
