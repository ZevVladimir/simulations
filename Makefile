all: nbody_sim hydro_sim adv_1d

nbody_sim:
	gcc -fopenmp -o nbody_sim.o nbody_sim.c nrutil.c integrators.c barnes_hut.c calc_fxns.c -lm

hydro_sim:
	gcc -fopenmp -o hydro_sim.o hydro_sim.c -lm

adv_1d:
	gcc -o adv_1d.o advection_1d.c -lm

clean:
	rm -f *.o nbody_sim hydro_sim