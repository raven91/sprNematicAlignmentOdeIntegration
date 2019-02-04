//
// Created by Nikita Kruk on 16.11.18.
//

#ifndef SPRNEMATICALIGNMENTODEINTEGRATION_BINARYOBSERVER_HPP
#define SPRNEMATICALIGNMENTODEINTEGRATION_BINARYOBSERVER_HPP

#include "../Definitions.hpp"
#include "../Parameters/PeriodicBoundaryConditionsConfiguration.hpp"
#include "../Parallelization/Thread.hpp"

#include <string>
#include <fstream>
#include <chrono>

class BinaryObserver
{
 public:

  explicit BinaryObserver(Thread &thread,
                          Real sigma,
                          Real rho,
                          Real D_T,
                          Real D_R,
                          PeriodicBoundaryConditionsConfiguration &pbc_config,
                          Real dt,
                          int trial);
  ~BinaryObserver();

  void SaveSystemState(std::vector<Real> &system_state, Real t);
  void SaveSummaryStatistics(std::vector<Real> &system_state, Real t);

 private:

  Thread &thread_;
  PeriodicBoundaryConditionsConfiguration &pbc_config_;
  std::chrono::time_point<std::chrono::system_clock> integration_step_timer_;
  int output_time_counter_[2];
  int output_time_threshold_[2];

  std::string simulation_file_name_;
  std::ofstream simulation_file_;

  std::string summary_statistics_file_name_;
  std::ofstream summary_statistics_file_;

};

#endif //SPRNEMATICALIGNMENTODEINTEGRATION_BINARYOBSERVER_HPP
