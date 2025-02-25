# HPQCGroupB
# How to run the ideal_gas file
-- locate it in the directory and run ./ideal_gas
-- If at any point you want to make changes to the code, you can compile the code using, gcc ideal_gas.c -o ideal_gas -lm
# How to run the ideal_gas_mpi file
-- locate it in the directory and run --mca btl tcp,self -np 4 ./ideal_gas_mpi
-- The amount of processes is controlled by the number 4
-- If at any point you want to make changes to the code, you can compile the code using, mpicc ideal_gas_mpi.c -o ideal_gas_mpi -lm
