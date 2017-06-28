CFLAGS = -Wall -std=gnu99 `gsl-config --cflags`
LDLIBS = $$(gsl-config --libs)
CC = gcc

all:	plot.png

plot.png: plot.gpi data
	gnuplot plot.gpi

data: main
	./main 1> data

main: main.o ode_driver.o rungekutta.o

.PHONEY: clean
clean:
	$(RM) main *.o *.txt data *.png
