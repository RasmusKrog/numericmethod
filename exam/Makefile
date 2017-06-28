CFLAGS = -Wall -std=gnu99 -O $$(gsl-config --cflags) -fopenmp
LDLIBS = $$(gsl-config --libs)
CC = gcc-4.9

all:	data pseudo.png quasi.png halton.png

pseudo.png: pseudo.gpi pseudo.dat
	gnuplot pseudo.gpi

quasi.png: quasi.gpi quasi.dat
	gnuplot quasi.gpi

halton.png: quasi.gpi halton.dat
	gnuplot quasi.gpi

data: main Makefile
	./main > data

pseudo.dat: points
		./points pseudo > pseudo.dat
quasi.dat: points
		./points quasi > quasi.dat
halton.dat: points
		./$< halton > $@

main: main.o lattice.o plainmc.o pseudo_random_vector.o halton.o quasi_random_vector.o
	$(CC) -o $@ $^ $(CFLAGS) $(LDLIBS)

points: points.o pseudo_random_vector.o quasi_random_vector.o halton.o 
	$(CC) -o $@ $^ $(CFLAGS) $(LDLIBS)

.PHONY: clean
clean: 
	$(RM) main *.o *.dat data *.png
