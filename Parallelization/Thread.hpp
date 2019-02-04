//
// Created by Nikita Kruk on 16.11.18.
//

#ifndef SPRNEMATICALIGNMENTODEINTEGRATION_THREAD_HPP
#define SPRNEMATICALIGNMENTODEINTEGRATION_THREAD_HPP

#include "../Definitions.hpp"

#include <vector>

class Thread
{
 public:

  Thread(int argc, char **argv);
  ~Thread();

  void SynchronizeVectorThroughBuffer(std::vector<Real> &vec);
  bool IsRoot();
  void BroadcastVector(std::vector<Real> &vec);
  const std::vector<int> &GetLoopIndices();
  void BroadcastCondition(int &condition);
  int GetNumberOfMpichThreads();
  int GetRank();

 private:

  int root_rank_;
  int rank_;
  int number_of_mpich_threads_;
  int number_of_elements_per_mpich_thread_;
  std::vector<int> loop_indices_;

};

#endif //SPRNEMATICALIGNMENTODEINTEGRATION_THREAD_HPP
