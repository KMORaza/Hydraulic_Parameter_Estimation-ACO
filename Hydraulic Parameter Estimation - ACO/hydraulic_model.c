#include "aco_hydraulic.h"
#include <math.h>
#define GRAVITY 9.81
#define KINEMATIC_VISCOSITY 1e-6
double calculate_head_loss(double flow, double diameter, double length,
                         double roughness, FlowEquationType eq_type) {
    if (diameter <= 0 || length <= 0) return 0.0;
    double head_loss = 0.0;
    double area = PI * pow(diameter/2, 2);
    if (area <= 0) return 0.0;
    double velocity = fabs(flow) / area;
    switch(eq_type) {
        case HAZEN_WILLIAMS:
            if (roughness > 0) {
                head_loss = 10.67 * pow(fabs(flow), 1.852) * length /
                          (pow(roughness, 1.852) * pow(diameter, 4.871));
            }
            break;
        case DARCY_WEISBACH: {
            double Re = velocity * diameter / KINEMATIC_VISCOSITY;
            double f;
            if (Re > 4000) {
                // Colebrook-White approximation
                double eps_over_D = roughness / diameter;
                f = 0.25 / pow(log10(eps_over_D/3.7 + 5.74/pow(Re,0.9)), 2);
            } else {
                f = 64.0 / Re;
            }
            head_loss = f * length * pow(velocity, 2) / (diameter * 2 * GRAVITY);
            break;
        }
        case MANNING:
            if (roughness > 0) {
                double Rh = diameter / 4.0;
                head_loss = pow(velocity * roughness / pow(Rh, 2.0/3.0), 2) * length;
            }
            break;
    }
    return head_loss;
}
void hydraulic_simulation(ACOConfig *config, const double params[]) {
    for (int i = 0; i < config->num_pipes; i++) {
        config->pipes[i].roughness = params[i];
    }
    const double source_head = 100.0;
    const double tolerance = 0.01;
    const int max_iter = 50;
    for (int i = 0; i < config->num_nodes; i++) {
        if (config->nodes[i].id == 1) {
            config->nodes[i].simulated_pressure = GRAVITY * (source_head - config->nodes[i].elevation);
        } else {
            config->nodes[i].simulated_pressure = GRAVITY * (source_head * 0.9 - config->nodes[i].elevation);
        }
    }
    for (int iter = 0; iter < max_iter; iter++) {
        double max_diff = 0.0;

        for (int i = 0; i < config->num_nodes; i++) {
            if (config->nodes[i].id == 1) continue;
            double sum_flow = 0.0;
            double sum_flow_head = 0.0;
            for (int j = 0; j < config->num_pipes; j++) {
                if (config->pipes[j].to_node == config->nodes[i].id) {
                    int upstream_idx = 0;
                    while (upstream_idx < config->num_nodes &&
                           config->nodes[upstream_idx].id != config->pipes[j].from_node) {
                        upstream_idx++;
                    }
                    if (upstream_idx < config->num_nodes) {
                        double upstream_head = config->nodes[upstream_idx].simulated_pressure/GRAVITY +
                                             config->nodes[upstream_idx].elevation;
                        double downstream_head = config->nodes[i].simulated_pressure/GRAVITY +
                                               config->nodes[i].elevation;
                        double delta_H = upstream_head - downstream_head;
                        double flow_sign = delta_H >= 0 ? 1.0 : -1.0;
                        double hl = calculate_head_loss(config->pipes[j].flow,
                                                       config->pipes[j].diameter,
                                                       config->pipes[j].length,
                                                       config->pipes[j].roughness,
                                                       config->flow_eq);
                        double flow_est = config->pipes[j].flow * sqrt(fabs(delta_H)/fmax(hl, 1e-6));
                        sum_flow += flow_est * flow_sign;
                        sum_flow_head += flow_est * flow_sign * upstream_head;
                    }
                }
            }
            sum_flow -= config->nodes[i].demand;
            if (fabs(sum_flow) > 1e-6) {
                double new_head = sum_flow_head / sum_flow;
                double new_pressure = GRAVITY * (new_head - config->nodes[i].elevation);
                double diff = fabs(new_pressure - config->nodes[i].simulated_pressure);
                if (diff > max_diff) max_diff = diff;
                config->nodes[i].simulated_pressure = new_pressure;
            }
        }
        if (max_diff < tolerance) break;
    }
}
double evaluate_fitness(ACOConfig *config) {
    double total_error = 0.0;
    int valid_nodes = 0;
    for (int i = 0; i < config->num_nodes; i++) {
        if (config->nodes[i].observed_pressure > 0 && config->nodes[i].id != 1) {
            double error = (config->nodes[i].simulated_pressure -
                          config->nodes[i].observed_pressure) /
                         config->nodes[i].observed_pressure;
            total_error += error * error;
            valid_nodes++;
        }
    }
    return valid_nodes > 0 ? sqrt(total_error / valid_nodes) * 100.0 : DBL_MAX;
}
