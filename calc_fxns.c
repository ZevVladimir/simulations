#include <math.h>
#include <stdio.h>

double calc_COM(double m1, double m2, double x1, double x2){
    // printf("%f,%f,%f,%f\n",m1,m2,x1,x2);
    return (m1 * x1 + m2 * x2) / (m1 + m2);
}

double calc_r(double x, double y){
    return sqrt(x * x + y * y);
}