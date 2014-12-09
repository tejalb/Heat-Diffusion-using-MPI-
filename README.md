Heat-Diffusion-using-MPI-
=========================

Solving the heat diffusion equation on a grid using OpenMP and MPI

Consider solving the heat diffusion equation with κ = constant = 1 on a two-dimensional domain of size 0 ≤ x, y ≤ π. Let the
boundary conditions be

T(x, 0) = cos^2 x
T(x, π) = sin^2 x
T(0, y) = T(π, y) (periodic in x)

This equation can be solved by centered finite differences in space and the
forward Euler method in time. This package includes 3 different
implementations:


• Serial: For the serial version, use heat_serial that runs
with command line options ./heat_serial {nx} for a solution with grid
size nx^2

• OpenMP: Parallel version heat_omp that
runs with command line options ./heat_omp {nx} {nthreads}. 

• MPI: Parallel version heat_mpi that runs
with mpiexec ./heat_mpi {nx}. Parallelized using domain decomposition.

