#ifndef ACO_HYDRAULIC_H
#define ACO_HYDRAULIC_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#define MAX_ANTS 50
#define MAX_PARAMS 12
#define MAX_NODES 50
#define MAX_PIPES 50
#define MAX_ITERATIONS 1000
#define MAX_STEPS 100

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
} ACOConfig;
void initialize_aco(ACOConfig *config, AntSolution ants[]);
void run_aco(ACOConfig *config, AntSolution ants[]);
double calculate_head_loss(double flow, double diameter, double length, double roughness, FlowEquationType eq_type);
double calculate_pressure(double head, double elevation);
void hydraulic_simulation(ACOConfig *config, double roughness_params[]);
double evaluate_fitness(ACOConfig *config);
double rand_range(double min, double max);
void init_random();
#endif
