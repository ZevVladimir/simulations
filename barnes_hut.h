#include "structs.h"
Node* init_node(double x_min, double x_max, double y_min, double y_max);
void subdivide(Node *node);
Node* build_tree(Particle *particles, int num_particles, double x_min, double x_max, double y_min, double y_max);
void free_tree(Node *node);
void calc_force(Node *node, Particle *particle, double theta, double epsilon, double G);
void clear_forces(Particle *particles, int num_ptls);
void print_tree(Node *node, int level, FILE *file);
void visualize_tree(Node *root, const char *filename);