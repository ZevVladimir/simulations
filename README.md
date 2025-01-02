# How to run
```
make
python3 gui.py
```
# How to use the GUI
Within the gui you can adjust the initial conditions of the galaxies as well as how many particles per galaxy, the number of snapshots, and the time step size. 

You then either:
- Click run simulation which will use the input parameters to run the nbody simulation and generate the data (depending on the number of particles and step number this can take a while). Once this is complete you can close the popup and run the video
- Click Play Existing Data which if you have already run the simulation before will just play that simulation back
- Click continue C simulation which if there is already data in the /data folder will continue the simulation from where the simulation had stopped
- Save to GIF if there is already data in the /data folder it will save it to a gif within that folder

At the top there is a slider which will dynamically adjust playback speed of the simulation.


# Files included

- .gitignore			ignore data files
- simulation.c         c code for the simulation 
- nbody.h 			Header file containing some stuff to aid computational work. 
- integrators.c 		Same integrators file we've been using in class. Using RK4. 
- integrators.h			Header for integrators.c
- nrutil.c 			
- nrutil.h 
- particle.h			Another header to help calculation work. 
- Makefile          compiles C code
- run_sim.sh        test script for running simulation without dealing with gui
- sim_gui.py            run this to interface with the code and get a nice video of the simulation
- particle.h        header file for particle containing constants
- BL_JA_ZM_ZV_writeup.pdf   The writeup for our project

