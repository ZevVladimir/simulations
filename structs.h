#ifndef PARTICLE_H
#define PARTICLE_H

// Structure to represent a single particle
typedef struct Particle {
    double x, y;     // Position in 2D
    double vx, vy;   // Velocity in 2D
    double fx, fy;   // Force in 2D
    double mass;     // Mass of the particle
} Particle;

typedef struct Node {
    // Bounds of the node's box
    double x_min, x_max;
    double y_min, y_max;

    // Total mass within the node and the location of the center of mass
    double tot_mass;
    double center_x, center_y;

    Particle *particle; // Pointer to particles within the node?
    struct Node *children[4]; // Child nodes
} Node;

// Optional: You can define any global constants related to particles here     // Gravitational constant (normalized)
static const double GRAVITATIONAL_CONSTANT = 1.0; // Softening parameter for gravity

#endif // PARTICLE_H

