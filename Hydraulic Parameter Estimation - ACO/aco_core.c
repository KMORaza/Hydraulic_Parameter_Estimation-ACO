#include "aco_hydraulic.h"
void initialize_aco(ACOConfig *config, AntSolution ants[]) {
    for (int i = 0; i < config->num_ants; i++) {
        for (int j = 0; j < config->num_params; j++) {
            double range = config->bounds[j].max - config->bounds[j].min;
            int steps = (int)(range / config->bounds[j].discretization);
            steps = steps > MAX_STEPS ? MAX_STEPS : steps;
            int random_step = rand() % (steps + 1);
            ants[i].params[j] = config->bounds[j].min + random_step * config->bounds[j].discretization;
        }
        ants[i].fitness = DBL_MAX;
        ants[i].probability = 0.0;
    }
}
void update_pheromones(double pheromone[MAX_PARAMS][MAX_STEPS], ACOConfig *config,
                      AntSolution ants[], double global_best_fitness,
                      const double global_best_params[]) {
    for (int j = 0; j < config->num_params; j++) {
        int steps = (int)((config->bounds[j].max - config->bounds[j].min) /
                         config->bounds[j].discretization) + 1;
        steps = steps > MAX_STEPS ? MAX_STEPS : steps;

        for (int k = 0; k < steps; k++) {
            pheromone[j][k] *= (1.0 - config->evaporation_rate);
        }
    }
    for (int i = 0; i < config->num_ants; i++) {
        for (int j = 0; j < config->num_params; j++) {
            int steps = (int)((config->bounds[j].max - config->bounds[j].min) /
                             config->bounds[j].discretization) + 1;
            steps = steps > MAX_STEPS ? MAX_STEPS : steps;

            int k = (int)((ants[i].params[j] - config->bounds[j].min) /
                         config->bounds[j].discretization + 0.5);
            k = k < 0 ? 0 : (k >= steps ? steps - 1 : k);
            pheromone[j][k] += 1.0 / ants[i].fitness;
        }
    }
    for (int j = 0; j < config->num_params; j++) {
        int steps = (int)((config->bounds[j].max - config->bounds[j].min) /
                         config->bounds[j].discretization) + 1;
        steps = steps > MAX_STEPS ? MAX_STEPS : steps;

        int k = (int)((global_best_params[j] - config->bounds[j].min) /
                     config->bounds[j].discretization + 0.5);
        k = k < 0 ? 0 : (k >= steps ? steps - 1 : k);
        pheromone[j][k] += 2.0 / global_best_fitness;
    }
}
void run_aco(ACOConfig *config, AntSolution ants[]) {
    double pheromone[MAX_PARAMS][MAX_STEPS] = {0};
    double heuristic[MAX_PARAMS][MAX_STEPS] = {0};
    double global_best_fitness = DBL_MAX;
    double global_best_params[MAX_PARAMS] = {0};
    for (int i = 0; i < config->num_params; i++) {
        int steps = (int)((config->bounds[i].max - config->bounds[i].min) /
                         config->bounds[i].discretization) + 1;
        steps = steps > MAX_STEPS ? MAX_STEPS : steps;
        for (int j = 0; j < steps; j++) {
            pheromone[i][j] = config->tau0;
            double param_value = config->bounds[i].min + j * config->bounds[i].discretization;
            double middle = (config->bounds[i].max + config->bounds[i].min) / 2.0;
            heuristic[i][j] = 1.0 / (1.0 + fabs(param_value - middle));
        }
    }
    for (int iter = 0; iter < config->max_iterations; iter++) {
        for (int i = 0; i < config->num_ants; i++) {
            hydraulic_simulation(config, ants[i].params);
            ants[i].fitness = evaluate_fitness(config);
            if (ants[i].fitness < global_best_fitness) {
                global_best_fitness = ants[i].fitness;
                memcpy(global_best_params, ants[i].params, sizeof(double) * config->num_params);
            }
        }
        update_pheromones(pheromone, config, ants, global_best_fitness, global_best_params);
        for (int i = 0; i < config->num_ants; i++) {
            for (int j = 0; j < config->num_params; j++) {
                int steps = (int)((config->bounds[j].max - config->bounds[j].min) /
                                 config->bounds[j].discretization) + 1;
                steps = steps > MAX_STEPS ? MAX_STEPS : steps;
                double total = 0.0;
                double probabilities[MAX_STEPS] = {0};
                for (int k = 0; k < steps; k++) {
                    probabilities[k] = pow(pheromone[j][k], config->alpha) *
                                     pow(heuristic[j][k], config->beta);
                    total += probabilities[k];
                }
                for (int k = 0; k < steps; k++) {
                    probabilities[k] /= total;
                }
                if (rand_range(0.0, 1.0) < config->q0) {
                    int best_k = 0;
                    double best_prob = 0.0;
                    for (int k = 0; k < steps; k++) {
                        if (probabilities[k] > best_prob) {
                            best_prob = probabilities[k];
                            best_k = k;
                        }
                    }
                    ants[i].params[j] = config->bounds[j].min + best_k * config->bounds[j].discretization;
                } else {
                    int selected_k = roulette_wheel_selection(probabilities, steps);
                    ants[i].params[j] = config->bounds[j].min + selected_k * config->bounds[j].discretization;
                }
            }
        }
        if (iter % 100 == 0) {
            printf("Iteration %4d | Best Fitness: %.6f\n", iter, global_best_fitness);
        }
    }
    printf("\n=== Optimization Results ===\n");
    printf("Best Fitness: %.6f\n", global_best_fitness);
    printf("Optimal Parameters:\n");
    for (int i = 0; i < config->num_params; i++) {
        printf("  Pipe %2d Roughness: %.4f\n", i+1, global_best_params[i]);
    }
    hydraulic_simulation(config, global_best_params);
    printf("\nNode Pressures (Simulated vs Observed):\n");
    for (int i = 0; i < config->num_nodes; i++) {
        printf("  Node %2d: %.2f kPa (sim) vs %.2f kPa (obs)\n",
               config->nodes[i].id,
               config->nodes[i].simulated_pressure,
               config->nodes[i].observed_pressure);
    }
}
