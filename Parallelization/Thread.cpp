//
// Created by Nikita Kruk on 16.11.18.
//

#include "Thread.hpp"

#include <cassert>
#include <iostream>
#include <numeric>  // std::iota

#if defined(MPI_PARALLELIZATION_IN_XY)
#include <mpi.h>
#endif

#if defined(MPI_PARALLELIZATION_IN_XY)
Thread::Thread(int argc, char **argv)
{
  root_rank_ = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &number_of_mpich_threads_);
  assert(!(kN % number_of_mpich_threads_));  // width by height must be divisible by the number of threads for this kind of parallelization
  number_of_elements_per_mpich_thread_ = kN / number_of_mpich_threads_;
}
#else
Thread::Thread(int argc, char **argv)
{
  root_rank_ = 0;
  rank_ = 0;
  number_of_mpich_threads_ = 1;
  assert(!(kN
      % number_of_mpich_threads_));  // width by height must be divisible by the number of threads for this kind of parallelization
  number_of_elements_per_mpich_thread_ = kN / number_of_mpich_threads_;

  loop_indices_ = std::vector<int>(number_of_elements_per_mpich_thread_, 0);
  std::iota(loop_indices_.begin(), loop_indices_.end(), rank_ * number_of_elements_per_mpich_thread_);
}
#endif

#if defined(MPI_PARALLELIZATION_IN_XY)
Thread::~Thread()
{

}
#else
Thread::~Thread()
{

}
#endif

#if defined(MPI_PARALLELIZATION_IN_XY)
void Thread::SynchronizeVectorThroughBuffer(std::vector<Real> &vec)
{
  std::vector<Real> buf(kN, 0.0);

  // TODO: link Real to MPI_DATA_TYPE
  MPI_Gather(&vec[rank_ * number_of_elements_per_mpich_thread_ * kS],
             number_of_elements_per_mpich_thread_ * kS,
             MPI_DOUBLE,
             &buf[0],
             number_of_elements_per_mpich_thread_ * kS,
             MPI_DOUBLE,
             root_rank_,
             MPI_COMM_WORLD);
  if (rank_ == root_rank_)
  {
    vec.swap(buf);
  }
  MPI_Bcast(&vec[0], (int) vec.size(), MPI_DOUBLE, root_rank_, MPI_COMM_WORLD);
}
#else
void Thread::SynchronizeVectorThroughBuffer(std::vector<Real> &vec)
{

}
#endif

#if defined(MPI_PARALLELIZATION_IN_XY)
bool Thread::IsRoot()
{
  return (rank_ == root_rank_);
}
#else
bool Thread::IsRoot()
{
  return true;
}
#endif

#if defined(MPI_PARALLELIZATION_IN_XY)
void Thread::BroadcastVector(std::vector<Real> &vec)
{
  MPI_Bcast(&vec[0], (int) vec.size(), MPI_DOUBLE, root_rank_, MPI_COMM_WORLD);
}
#else
void Thread::BroadcastVector(std::vector<Real> &vec)
{

}
#endif

#if defined(MPI_PARALLELIZATION_IN_XY)
const std::vector<int> &Thread::GetLoopIndices()
{
  return loop_indices_;
}
#else
const std::vector<int> &Thread::GetLoopIndices()
{
  return loop_indices_;
}
#endif

#if defined(MPI_PARALLELIZATION_IN_XY)
void Thread::BroadcastCondition(int &condition)
{
  MPI_Bcast(&condition, 1, MPI_INT, root_rank_, MPI_COMM_WORLD);
}
#else
void Thread::BroadcastCondition(int &condition)
{

}
#endif

#if defined(MPI_PARALLELIZATION_IN_XY)
int Thread::GetNumberOfMpichThreads()
{
  return number_of_mpich_threads_;
}
#else
int Thread::GetNumberOfMpichThreads()
{
  return number_of_mpich_threads_;
}
#endif

#if defined(MPI_PARALLELIZATION_IN_XY)
int Thread::GetRank()
{
  return rank_;
}
#else
int Thread::GetRank()
{
  return rank_;
}
#endif