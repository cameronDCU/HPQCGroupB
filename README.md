# HPQCGroupB
# How to run the ideal_gas file
-- locate it in the directory and run ./ideal_gas <num_particles>
-- If at any point you want to make changes to the code, you can open the code using, nano ideal_gas.c, and compile the code using, gcc ideal_gas.c -o ideal_gas -lm
# How to run the ideal_gas_mpi file
-- locate it in the directory and run --mca btl tcp,self -np <num_processes> ./ideal_gas_mpi <num_processes>
-- If at any point you want to make changes to the code, you can open the code using nano ideal_gas_mpi.c, and compile the code using, mpicc ideal_gas_mpi.c -o ideal_gas_mpi -lm
