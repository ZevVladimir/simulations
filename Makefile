all: nbody_sim hydro_sim adv_1d adv_2d

nbody_sim:
	gcc -fopenmp -o nbody_sim.o nbody_sim.c nrutil.c integrators.c barnes_hut.c calc_fxns.c -lm

hydro_sim:
	gcc -fopenmp -o hydro_sim.o hydro_sim.c -lm

adv_1d:
	gcc -o adv_1d.o advection_1d.c -lm

adv_2d:
	gcc -o adv_2d.o advection_2d.c -lm

room_2d:
	gcc -I. -o room_2d room_2d.c Navier_Stokes/linalg.c Navier_Stokes/fdm.c -lm

clean:
	rm -f *.o nbody_sim hydro_sim