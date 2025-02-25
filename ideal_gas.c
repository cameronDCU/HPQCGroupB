#include <stdio.h>   // For printf
#include <stdlib.h>  // For rand, srand, malloc, free
#include <math.h>    // For sqrt, log, cos
#include <time.h>    // For clock(), CLOCKS_PER_SEC

#define L 10.0  // Box size
#define T 1.0   // Temperature scale for velocity distribution
#define DT 0.01 // Time step
#define STEPS 1000  // Number of simulation steps
#define MASS 1.0 // Particle mass

typedef struct {
    double x, y, z;
} Vector3;

// Generate a uniform random number between [0, 1]
double rand_uniform() {
    return rand() / (double)RAND_MAX;
}

// Generate a normally distributed random number using Box-Muller transform
double rand_normal() {
    double u1 = rand_uniform();
    double u2 = rand_uniform();
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

// Initialize positions uniformly in the box
void initialize_positions(Vector3 *pos, int num_particles) {
    for (int i = 0; i < num_particles; i++) {
        pos[i].x = rand_uniform() * L;
        pos[i].y = rand_uniform() * L;
        pos[i].z = rand_uniform() * L;
    }
}

// Initialize velocities using Maxwell-Boltzmann distribution
void initialize_velocities(Vector3 *vel, int num_particles) {
    for (int i = 0; i < num_particles; i++) {
        double v_mag = sqrt(-2.0 * T * log(rand_uniform()));
        Vector3 v_dir = {rand_normal(), rand_normal(), rand_normal()};
        double norm = sqrt(v_dir.x * v_dir.x + v_dir.y * v_dir.y + v_dir.z * v_dir.z);
        if (norm == 0) norm = 1.0;
        vel[i].x = (v_mag * v_dir.x) / norm;
        vel[i].y = (v_mag * v_dir.y) / norm;
        vel[i].z = (v_mag * v_dir.z) / norm;
    }
}

// Update positions and handle wall collisions
void update_positions(Vector3 *pos, Vector3 *vel, int num_particles) {
    for (int i = 0; i < num_particles; i++) {
        pos[i].x += vel[i].x * DT;
        pos[i].y += vel[i].y * DT;
        pos[i].z += vel[i].z * DT;

        // Wall collision handling
        if (pos[i].x < 0 || pos[i].x > L) {
            vel[i].x *= -1;
            pos[i].x = (pos[i].x < 0) ? 0 : L;
        }
        if (pos[i].y < 0 || pos[i].y > L) {
            vel[i].y *= -1;
            pos[i].y = (pos[i].y < 0) ? 0 : L;
        }
        if (pos[i].z < 0 || pos[i].z > L) {
            vel[i].z *= -1;
            pos[i].z = (pos[i].z < 0) ? 0 : L;
        }
    }
}

// Compute temperature, net velocity, and pressure
void compute_properties(Vector3 *vel, int num_particles, double *temperature, double *net_velocity, double *pressure) {
    double kinetic_energy = 0.0;
    Vector3 momentum = {0.0, 0.0, 0.0};
    
    for (int i = 0; i < num_particles; i++) {
        double v_sq = vel[i].x * vel[i].x + vel[i].y * vel[i].y + vel[i].z * vel[i].z;
        kinetic_energy += 0.5 * MASS * v_sq;
        momentum.x += vel[i].x;
        momentum.y += vel[i].y;
        momentum.z += vel[i].z;
    }

    *temperature = (2.0 / 3.0) * (kinetic_energy / num_particles);
    *net_velocity = sqrt(momentum.x * momentum.x + momentum.y * momentum.y + momentum.z * momentum.z) / num_particles;
    *pressure = (2.0 * kinetic_energy) / (3.0 * L * L * L);
}

int main(int argc, char *argv[]) {
    clock_t start_time, end_time;
    start_time = clock();  // Start timing execution

    srand(time(NULL)); // Seed random number generator

    if (argc != 2) {fprintf(stderr, "Usage: %s <number_of_particles>\n", argv[0]);
        return 1;
    }

       // Convert input argument to integer
    int N = atoi(argv[1]);
    if (N <= 0) {
        fprintf(stderr, "Error: Number of particles must be a positive integer.\n");
        return 1;
    }


    Vector3 *pos = (Vector3 *)malloc(N * sizeof(Vector3));
    Vector3 *vel = (Vector3 *)malloc(N * sizeof(Vector3));

    initialize_positions(pos, N);
    initialize_velocities(vel, N);

    for (int step = 0; step < STEPS; step++) {
        update_positions(pos, vel, N);
    }

    double temperature, net_velocity, pressure;
    compute_properties(vel, N, &temperature, &net_velocity, &pressure);

    printf("System Properties:\n");
    printf("Temperature: %.2f\n", temperature);
    printf("Net Velocity: %.2f\n", net_velocity);
    printf("Pressure: %.2f\n", pressure);

    free(pos);
    free(vel);

    end_time = clock();  // Stop timing execution
    double execution_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Execution Time: %.6f seconds\n", execution_time);

    return 0;
}
