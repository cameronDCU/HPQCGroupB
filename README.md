# HPQCGroupB
## How to run the programs
 
### How to run the ideal_gas file
-- locate it in the directory and run ./ideal_gas <num_particles>
-- If at any point you want to make changes to the code, you can open the code using, nano ideal_gas.c, and compile the code using, gcc ideal_gas.c -o ideal_gas -lm
### How to run the ideal_gas_mpi file
-- locate it in the directory and run --mca btl tcp,self -np <num_processes> ./ideal_gas_mpi <num_processes>
-- If at any point you want to make changes to the code, you can open the code using nano ideal_gas_mpi.c, and compile the code using, mpicc ideal_gas_mpi.c -o ideal_gas_mpi -lm
### How to run the ideal_gas.py program
-- Run the line "python ideal_gas.py <num_processes>"

## How The Programs Work

The simulation initializes particles with random positions and velocities following the Maxwell-Boltzmann distribution, updates their motion over time, and handles elastic wall collisions. Key macroscopic properties such as temperature, pressure, and net velocity are computed to verify system behavior. The Python implementation uses NumPy for efficiency, while the C version leverages direct memory management for improved performance. The MPI version distributes particles across processes, updating them independently and aggregating results using MPI_Reduce, enabling scalability for large simulations. The project highlights a transition from an easy-to-develop Python model to a highly efficient parallelized approach in C.


## Benchmark

The codes were ran for various timesteps and the results can be seen below.

![Benchmarking](https://github.com/cameronDCU/HPQCGroupB/blob/main/benchmark_group_hpqc.png)

As can be observed, the C implementations of the code outperformed the Python prototype 'ideal_gas.py'. As the number of particles increased, MPI became more suitable for the task over the serial C implementation.
