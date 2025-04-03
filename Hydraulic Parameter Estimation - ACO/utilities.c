#include "aco_hydraulic.h"
void init_random() {
    srand((unsigned int)time(NULL));
}
double rand_range(double min, double max) {
    return min + (max - min) * ((double)rand() / RAND_MAX);
}
int roulette_wheel_selection(const double probabilities[], int n) {
    double r = rand_range(0.0, 1.0);
    double sum = 0.0;

    for (int i = 0; i < n; i++) {
        sum += probabilities[i];
        if (r <= sum) {
            return i;
        }
    }
    return n - 1;
}
void print_network_summary(ACOConfig *config) {
    printf("\n=== Network Summary ===\n");
    printf("Nodes: %d\n", config->num_nodes);
    printf("Pipes: %d\n", config->num_pipes);
    printf("Flow Equation: %s\n",
           config->flow_eq == HAZEN_WILLIAMS ? "Hazen-Williams" :
           config->flow_eq == DARCY_WEISBACH ? "Darcy-Weisbach" : "Manning");
    printf("\nNodes:\n");
    for (int i = 0; i < config->num_nodes; i++) {
        printf("  Node %2d: Elev=%.2fm, Demand=%.3fm3/s\n",
               config->nodes[i].id,
               config->nodes[i].elevation,
               config->nodes[i].demand);
    }
    printf("\nPipes:\n");
    for (int i = 0; i < config->num_pipes; i++) {
        printf("  Pipe %2d: %2d->%2d, L=%.1fm, D=%.3fm, Q=%.3fm3/s\n",
               config->pipes[i].id,
               config->pipes[i].from_node,
               config->pipes[i].to_node,
               config->pipes[i].length,
               config->pipes[i].diameter,
               config->pipes[i].flow);
    }
}
