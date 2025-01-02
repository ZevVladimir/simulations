void rk4(double y[], double dydx[], int n, double x, double h, double yout[], int N, int glxy1_N,
        void (*derivs)(double, double [], double [], int, int));
        
void lf2(double y[],double dydt[],int n,int N,double t,double h, double mass[],
        void (*derivs)(double [],double [],int,double []));

void lf2_BH(Node *node, Particle *particles, int n_ptl, double theta, double epsilon, double G, double h);