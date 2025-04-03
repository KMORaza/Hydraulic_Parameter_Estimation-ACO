#ifndef ACO_HYDRAULIC_H
#define ACO_HYDRAULIC_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <string.h>
#define MAX_ANTS 50
#define MAX_PARAMS 20
#define MAX_NODES 100
#define MAX_PIPES 100
#define MAX_ITERATIONS 1000
#define MAX_STEPS 100
#define PI 3.141592653589793

typedef enum {
    HAZEN_WILLIAMS,
    DARCY_WEISBACH,
    MANNING
} FlowEquationType;
typedef struct {
    int id;
    double elevation;
    double demand;
    double observed_pressure;
    double simulated_pressure;
} Node;
typedef struct {
    int id;
    int from_node;
    int to_node;
    double length;
    double diameter;
    double flow;
    double roughness;
} Pipe;
typedef struct {
    double min;
    double max;
    double discretization;
} ParamBounds;
typedef struct {
    double params[MAX_PARAMS];
    double fitness;
    double probability;
} AntSolution;
typedef struct {
    int num_ants;
    int num_params;
    int max_iterations;
    double evaporation_rate;
    double alpha;
    double beta;
    double q0;
    double tau0;
    FlowEquationType flow_eq;
    ParamBounds bounds[MAX_PARAMS];
    Node nodes[MAX_NODES];
    Pipe pipes[MAX_PIPES];
    int num_nodes;
    int num_pipes;
    double target_pressure;
} ACOConfig;
void initialize_aco(ACOConfig *config, AntSolution ants[]);
void run_aco(ACOConfig *config, AntSolution ants[]);
void update_pheromones(double pheromone[MAX_PARAMS][MAX_STEPS], ACOConfig *config, AntSolution ants[], double global_best_fitness, const double global_best_params[]);
double calculate_head_loss(double flow, double diameter, double length, double roughness, FlowEquationType eq_type);
void hydraulic_simulation(ACOConfig *config, const double roughness_params[]);
double evaluate_fitness(ACOConfig *config);
double rand_range(double min, double max);
int roulette_wheel_selection(const double probabilities[], int n);
void init_random();
void print_network_summary(ACOConfig *config);
#endif
