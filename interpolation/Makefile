CFLAGS = -Wall -std=gnu99 -O `gsl-config --cflags` 
LDLIBS = $$(gsl-config --libs)
CC = clang 

all:  lplot.png qplot.png cplot.png

lplot.png: lplot.gpi data
	gnuplot lplot.gpi

qplot.png: qplot.gpi data
	gnuplot qplot.gpi

cplot.png: cplot.gpi data
	gnuplot cplot.gpi

data: main Makefile
	./main 30 > data

main: main.o spline.o

.PHONEY: clean
clean: 
	$(RM) main *.o *.pdf *.txt data *.png
