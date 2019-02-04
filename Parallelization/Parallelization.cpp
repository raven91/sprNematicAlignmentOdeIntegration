//
// Created by Nikita Kruk on 16.11.18.
//

#include "Parallelization.hpp"

#if defined(MPI_PARALLELIZATION_IN_XY)
#include <mpi.h>
#endif

#if defined(MPI_PARALLELIZATION_IN_XY)
void LaunchParallelSession(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
}
#else
void LaunchParallelSession(int argc, char **argv)
{

}
#endif

#if defined(MPI_PARALLELIZATION_IN_XY)
void FinalizeParallelSession()
{
  MPI_Finalize();
}
#else
void FinalizeParallelSession()
{

}
#endif