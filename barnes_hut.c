#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structs.h"       // For Particle structure and constants
#include "integrators.h"    // For rk4 integrator
#include "calc_fxns.h"


// Instantiate a node
Node* init_node(double x_min, double x_max, double y_min, double y_max){
    Node *node = (Node*)malloc(sizeof(Node));
    node->x_min = x_min;
    node->x_max = x_max;
    node->y_min = y_min;
    node->y_max = y_max;
    node->tot_mass = 0.0;
    node->center_x = 0.0;
    node->center_y = 0.0;
    node->particle = NULL;
    for (int i = 0; i < 4; i++){
        node->children[i] = NULL;
    }
    return node;
}

// Split one node into 4 children nodes evenly
void subdivide(Node *node){
    double mid_x = (node->x_min + node->x_max) / 2.0; 
    double mid_y = (node->y_min + node->y_max) / 2.0; 
    node->children[0] = init_node(node->x_min,mid_x,node->y_min,mid_y); //BL
    node->children[1] = init_node(mid_x,node->x_max,node->y_min,mid_y); //BR
    node->children[2] = init_node(node->x_min,mid_x,mid_y,node->y_max); //TL
    node->children[3] = init_node(mid_x,node->x_max,mid_y,node->y_max); //TR
}

void insert_ptl(Node *node, Particle *particle) {
    // If the node is empty and has no child nodes
    if (node->particle == NULL && node->children[0] == NULL) {
        node->particle = particle;
        node->tot_mass = particle->mass;
        node->center_x = particle->x;
        node->center_y = particle->y;
        return;
    }

    // If the node isn't empty but does not have subdivisions
    if (node->children[0] == NULL) {
        subdivide(node);
    }

    // If the node contains a particle, move it into a child node
    if (node->particle != NULL) {
        Particle *existing = node->particle;
        node->particle = NULL;
        insert_ptl(node, existing);
    }

    // Determine which quadrant the particle belongs to
    double avg_x = (node->x_min + node->x_max) / 2.0;
    double avg_y = (node->y_min + node->y_max) / 2.0;
    int index = (particle->x > avg_x) + 2 * (particle->y > avg_y);
    insert_ptl(node->children[index], particle);

    // Recalculate total mass and center of mass
    node->tot_mass = 0.0;
    node->center_x = 0.0;
    node->center_y = 0.0;

    for (int i = 0; i < 4; i++) {
        if (node->children[i] != NULL) {
            Node *child = node->children[i];
            node->tot_mass += child->tot_mass;
            node->center_x += child->center_x * child->tot_mass;
            node->center_y += child->center_y * child->tot_mass;
        }
    }

    if (node->tot_mass > 0) {
        node->center_x /= node->tot_mass;
        node->center_y /= node->tot_mass;
    }
}


Node* build_tree(Particle *particles, int num_particles, double x_min, double x_max, double y_min, double y_max){
    Node *root = init_node(x_min, x_max, y_min, y_max); //create the root node
    // Then let insert_ptl take care of placing all the particles as subdivision is taken care of automatically
    for (int i = 0; i < num_particles; i++){
        insert_ptl(root, &particles[i]);
    }
    return root;
}

void free_tree(Node *node){
    // If the given node is empty reached the end of recursion
    if (node == NULL){
        return;
    }
    // Otherwise loop through each of the node's children recursively until the bottom one
    for (int i = 0; i < 4; i++){
        free_tree(node->children[i]);
    }
    // Free the node as recursion collapses
    free(node);
}

void calc_force(Node *node, Particle *particle, double theta, double epsilon, double G){
    if (node == NULL || node->tot_mass == 0){
        return;
    }

    // Calculate the node's size
    double size = node->x_max - node->x_min;

    // Calculate the distance between the particle and the node
    double dx = node->center_x - particle->x;
    double dy = node->center_y - particle->y;
    double r = calc_r(dx, dy);

    if (r == 0) {
        return; // Skip force calculation if the distance is zero
    }
    
    if((size / r < theta) || node->children[0] == NULL){
        double denom = r * r + epsilon * epsilon;

        if (denom == 0 || isnan(denom)) {
            return;
        }

        double F = (G * particle->mass * node->tot_mass) / (denom + sqrt(denom));

        particle->fx += F * (dx / r);
        particle->fy += F * (dy / r);
    }
    else{
        for(int i = 0; i < 4; i++){
            calc_force(node->children[i],particle,theta,epsilon,G);
        }
    }
}

void clear_forces(Particle *particles, int num_ptls){
    for (int i = 0; i < num_ptls; i++){
        particles[i].fx = 0.0;
        particles[i].fy = 0.0;
    }
}

void print_tree(Node *node, int level, FILE *file) {
    if (node == NULL) return;

    // Indentation for the current level
    for (int i = 0; i < level; i++) {
        fprintf(file, "  ");
    }

    // Print node details
    fprintf(file, "Node: Center of Mass (%.2f, %.2f), Total Mass: %.2f, X lim (%.3f, %.3f), Y lim (%.3f, %.3f)\n",
            node->center_x, node->center_y, node->tot_mass, node->x_min, node->x_max, node->y_min, node->y_max);

    // Recursively print children
    for (int i = 0; i < 4; i++) {
        if (node->children[i] != NULL) {
            print_tree(node->children[i], level + 1, file);
        }
    }
}

// Function to start the tree visualization
void visualize_tree(Node *root, const char *filename) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        printf("Error: Could not open file %s for writing.\n", filename);
        return;
    }

    print_tree(root, 0, file);
    fclose(file);
    printf("Tree visualization written to %s\n", filename);
}

