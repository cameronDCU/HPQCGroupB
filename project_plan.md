# Project Plan: Ideal Gas Simulation with MPI

## Objective
Convert the existing Python code for simulating an ideal gas in a 3D box to a C implementation using MPI. The particles should be distributed among 4 processes.

## Plan

**MPI Initialization**:
  - Initialize MPI with `MPI_Init()`.
  - Obtain the process rank using `MPI_Comm_rank()` and total number of processes using `MPI_Comm_size()`.
  - Set global parameters (box size, number of particles, simulation time step, etc.).

Use MPI_Scatter to distribute the particles among the processes.

Run simulation loop to update particle properties and check if a collision with the wall has occured.

Each process will calculate Temperature, Net velocity and Pressure.

Each process will send the results to the root process and sum them using MPI_Reduce.

The root process will return the properties of the system.

We will benchmark the performance of this program and compare it to a python program, and C program which do not use MPI.
