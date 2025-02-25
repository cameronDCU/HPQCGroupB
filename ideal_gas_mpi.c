#include <stdio.h>   // For printf
#include <stdlib.h>  // For rand, srand
#include <math.h>    // For sqrt, log, cos
#include <time.h>    // For time
#include <mpi.h>     // For MPI functions

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
            if (pos[i].x < 0) pos[i].x = 0;
            if (pos[i].x > L) pos[i].x = L;
        }
        if (pos[i].y < 0 || pos[i].y > L) {
            vel[i].y *= -1;
            if (pos[i].y < 0) pos[i].y = 0;
            if (pos[i].y > L) pos[i].y = L;
        }
        if (pos[i].z < 0 || pos[i].z > L) {
            vel[i].z *= -1;
            if (pos[i].z < 0) pos[i].z = 0;
            if (pos[i].z > L) pos[i].z = L;
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

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    srand(time(NULL) + rank); // Seed random number generator
    
    // Convert input argument to integer
    int N = atoi(argv[1]);
    if (N <= 0) {
        fprintf(stderr, "Error: Number of particles must be a positive integer.\n");
        return 1;
    }

    int local_N = N / size;
    Vector3 *local_pos = (Vector3 *)malloc(local_N * sizeof(Vector3));
    Vector3 *local_vel = (Vector3 *)malloc(local_N * sizeof(Vector3));

    // Start execution time measurement
    double start_time = MPI_Wtime();

    // Scatter positions
    Vector3 *global_pos = NULL;
    if (rank == 0) {
        global_pos = (Vector3 *)malloc(N * sizeof(Vector3));
        initialize_positions(global_pos, N);
    }

    MPI_Scatter(global_pos, local_N * sizeof(Vector3), MPI_BYTE, 
                local_pos, local_N * sizeof(Vector3), MPI_BYTE, 
                0, MPI_COMM_WORLD);

    if (rank == 0) {
        free(global_pos);
    }

    initialize_velocities(local_vel, local_N);

    for (int step = 0; step < STEPS; step++) {
        update_positions(local_pos, local_vel, local_N);
    }

    double local_temp, local_net_vel, local_pressure;
    compute_properties(local_vel, local_N, &local_temp, &local_net_vel, &local_pressure);

    double global_temp, global_net_vel, global_pressure;
    MPI_Reduce(&local_temp, &global_temp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_net_vel, &global_net_vel, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_pressure, &global_pressure, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Write the results to a CSV file (only at root process)
    if (rank == 0) {
        FILE *file = fopen("particle_properties.csv", "w");
        if (file == NULL) {
            fprintf(stderr, "Error opening file for writing\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        fprintf(file, "Temperature, Net_Velocity, Pressure\n");
        fprintf(file, "%.2f, %.2f, %.2f\n", global_temp, global_net_vel, global_pressure);
        fclose(file);
    }
// Stop execution time measurement
    double end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;

    if (rank == 0) {
        printf("System Properties:\n");
        printf("Temperature: %.2f\n", global_temp / size);
        printf("Net Velocity: %.2f\n", global_net_vel / size);
        printf("Pressure: %.2f\n", global_pressure / size);
        printf("Execution Time: %.6f seconds\n", elapsed_time);
    }

    free(local_pos);
    free(local_vel);
    MPI_Finalize();
    return 0;
}
