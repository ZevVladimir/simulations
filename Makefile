all: sim

sim:
	gcc -fopenmp -o sim.o simulation.c nrutil.c integrators.c barnes_hut.c calc_fxns.c -lm

clean:
	rm -f *.o sim