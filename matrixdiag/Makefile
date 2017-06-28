CFLAGS = -Wall -std=gnu99 -O `gsl-config --cflags`  #runs in #backticks,gives output back
LDLIBS = $$(gsl-config --libs)
CC = clang

all:  data
	cat data
	#plot.png

#plot.png: plot.gpi data
#	gnuplot plot.gpi

data: main Makefile
	./main > data

main: main.o gsl_shortcut.o matrixdiag.o

.PHONEY: clean
clean: 
	$(RM) main *.o *.txt data *.png
