Hydraulic Parameter Estimation using Ant Colony Optimization / Schätzung hydraulischer Parameter mittels Ameisenkolonieoptimierung
---

* Pipe roughness coefficients are represented as discrete values (nodes in a graph)
* Each ant's path represents a potential solution (set of roughness values)
* Initialization:
  * Define search bounds for each parameter
  * Create multiple "ants" with random initial solutions
  * Initialize pheromone trails uniformly
* Each ant selects roughness values probabilistically:
  * Favors values with higher pheromone (better historical performance)
  * Considers heuristic information (closeness to midpoint of range)
* Uses exploration-exploitation balance (parameter q₀)
* Hydraulic Evaluation:
  * Run hydraulic simulation with proposed roughness values
  * Calculate nodal pressures
  * Compare with observed pressures
  * Compute fitness (RMS % pressure error)
* Pheromone Update:
  * Evaporate all pheromone trails
  * Best solution gets extra pheromone
* Stop after maximum iterations or when solutions stabilize & return best-found roughness values
