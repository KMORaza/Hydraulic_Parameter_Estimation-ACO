#include "aco_hydraulic.h"
#include <stdio.h>
void run_optimization(FlowEquationType eq_type) {
    ACOConfig config = {
        .num_ants = 30,
        .num_params = 5,
        .max_iterations = 200,
        .evaporation_rate = 0.2,
        .alpha = 1.5,
        .beta = 2.0,
        .q0 = 0.9,
        .tau0 = 0.5,
        .flow_eq = eq_type,
        .bounds = {
            {120.0, 140.0, 1.0},
            {0.0001, 0.01, 0.0001},
            {0.01, 0.03, 0.001}
        }
    };
    if (eq_type == HAZEN_WILLIAMS) {
        for (int i = 0; i < config.num_params; i++) {
            config.bounds[i] = (ParamBounds){120.0, 140.0, 1.0};
        }
    } else if (eq_type == DARCY_WEISBACH) {
        for (int i = 0; i < config.num_params; i++) {
            config.bounds[i] = (ParamBounds){0.0001, 0.01, 0.0001};
        }
    } else {
        for (int i = 0; i < config.num_params; i++) {
            config.bounds[i] = (ParamBounds){0.01, 0.03, 0.001};
        }
    }
    config.num_nodes = 4;
    config.nodes[0] = (Node){1, 50.0, 0.0, 490.5};
    config.nodes[1] = (Node){2, 45.0, 0.15, 350.0};
    config.nodes[2] = (Node){3, 40.0, 0.20, 320.0};
    config.nodes[3] = (Node){4, 42.0, 0.10, 340.0};
    config.num_pipes = 5;
    config.pipes[0] = (Pipe){1, 1, 2, 100.0, 0.25, 0.30, 0.0};
    config.pipes[1] = (Pipe){2, 2, 3, 150.0, 0.20, 0.20, 0.0};
    config.pipes[2] = (Pipe){3, 2, 4, 120.0, 0.18, 0.10, 0.0};
    config.pipes[3] = (Pipe){4, 3, 4, 80.0, 0.15, 0.05, 0.0};
    config.pipes[4] = (Pipe){5, 1, 3, 200.0, 0.30, 0.10, 0.0};
    printf("\n=== Optimizing with %s ===",
           eq_type == HAZEN_WILLIAMS ? "Hazen-Williams" :
           eq_type == DARCY_WEISBACH ? "Darcy-Weisbach" : "Manning");
    print_network_summary(&config);
    AntSolution ants[MAX_ANTS];
    initialize_aco(&config, ants);
    run_aco(&config, ants);
}
int main() {
    init_random();
    run_optimization(HAZEN_WILLIAMS);
    run_optimization(DARCY_WEISBACH);
    run_optimization(MANNING);
    return 0;
}
