#define NRANSI
#include "nrutil.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "barnes_hut.h"
#include <omp.h>

void rk4(double y[], double dydx[], int n, double x, double h, double yout[], int N, int glxy1_N,
	void (*derivs)(double, double [], double [], int, int))
{
	int i;
	double xh,hh,h6;

	double *dym = (double *)malloc(n * sizeof(double)); 
	double *dyt = (double *)malloc(n * sizeof(double)); 
	double *yt = (double *)malloc(n * sizeof(double)); 
	// dym=dvector(1,n);
	// dyt=dvector(1,n);
	// yt=dvector(1,n);

	hh=h*0.5;
	h6=h/6.0;
	xh=x+hh;

	for (i=0;i<n;i++) yt[i]=y[i]+hh*dydx[i];
	(*derivs)(xh,yt,dyt,N,glxy1_N);
	for (i=0;i<n;i++) yt[i]=y[i]+hh*dyt[i];
	(*derivs)(xh,yt,dym,N,glxy1_N);
	for (i=0;i<n;i++) {
		yt[i]=y[i]+h*dym[i];
		dym[i] += dyt[i];
	}
	(*derivs)(x+h,yt,dyt,N,glxy1_N);
	for (i=0;i<n;i++)
		yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);

	// free_dvector(yt,1,n);
	// free_dvector(dyt,1,n);
	// free_dvector(dym,1,n);
	free(dym);
	free(dyt);
	free(yt);
}

void lf2(double y[],double dydt[],int n,int N,double t,double h, double mass[],
void (*derivs)(double [],double [],int,double []))
{
/*
** Performs one LeapFrog step (size h) on the n variables in y, using
** the SUPPLIED derivatives dydt at time t (the user must precompute
** the derivatives before the first call and should not modify dydt
** between calls). Also supply the function for computing the derivs,
** which will be called for the closing kick. The generic pointer p
** will be passed to the derivatives function. Note that the values
** in y will be time synchronized on return. Finally, note that this
** function is only appropriate for certain special classes of ODE.
*/
// printf("\n%f,%f,%f,%f,%f,%f\n",y[7],y[9],y[11],y[8],y[10],y[12]);
// printf("%f,%f,%f,%f,%f,%f\n",dydt[7],dydt[9],dydt[11],dydt[8],dydt[10],dydt[12]);
int i;

assert(n%2 == 0); /* must be even number of equations for leapfrog */
/* opening kick: "velocities" are in elements 2,4,...,n-1 */
// v_{n+1/2} =  v_n + h/2 * a
for (i=1;i<=n;i+=2){
	y[i+1] += 0.5*h*dydt[i+1];
}

/* drift: "positions" are in elements 1,3,...,n-2 */
// y_{n+1} = y_{n} + h * v_{n+1/2}
for (i=1;i<=n;i+=2){
	y[i] += h*y[i+1];
}

/* get new derivatives */
derivs(y,dydt,N,mass); /* dydt[1],dydt[3]....,dydt[n-2] not used */

/* closing kick */
for (i=1;i<=n;i+=2){
	y[i+1] += 0.5*h*dydt[i+1];
}

// printf("%f,%f,%f,%f,%f,%f\n",y[7],y[9],y[11],y[8],y[10],y[12]);
}


void lf2_BH(Node *node, Particle *particles, int n_ptl, double theta, double epsilon, double G, double h){
	clear_forces(particles,n_ptl);

	#pragma omp parallel for
	for (int i = 0; i < n_ptl; i++){
		// Calculate forces on this particle
		calc_force(node,&particles[i],theta,epsilon,G);

		// Update velocities 1/2 step
		particles[i].vx += 0.5 * particles[i].fx / particles[i].mass * h;
		particles[i].vy += 0.5 * particles[i].fy / particles[i].mass * h;

		// Update positions
		particles[i].x += particles[i].vx * h;
		particles[i].y += particles[i].vy * h;

		// Update velocities (second half-step)
		particles[i].vx += 0.5 * particles[i].fx / particles[i].mass * h;
		particles[i].vy += 0.5 * particles[i].fy / particles[i].mass * h;
	}
}

#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */
