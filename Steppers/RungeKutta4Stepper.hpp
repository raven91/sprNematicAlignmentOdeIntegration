//
// Created by Nikita Kruk on 16.11.18.
//

#ifndef SPRNEMATICALIGNMENTODEINTEGRATION_RUNGEKUTTA4STEPPER_HPP
#define SPRNEMATICALIGNMENTODEINTEGRATION_RUNGEKUTTA4STEPPER_HPP

#include "../Definitions.hpp"
#include "../DynamicalSystems/NematicAlignmentSystem.hpp"
#include "../DynamicalSystems/NormalizedVelocitySystem.hpp"
#include "../Parallelization/Thread.hpp"

class RungeKutta4Stepper
{
 public:

  explicit RungeKutta4Stepper(Thread &thread);
  ~RungeKutta4Stepper();

  void DoStep(NematicAlignmentSystem &system, std::vector<Real> &system_state, Real t, Real dt);
  void DoStep(NormalizedVelocitySystem &system, std::vector<Real> &system_state, Real t, Real dt);

 private:

  Thread &thread_;
  std::vector<Real> k_1_;
  std::vector<Real> k_2_;
  std::vector<Real> k_3_;
  std::vector<Real> k_4_;

};

#endif //SPRNEMATICALIGNMENTODEINTEGRATION_RUNGEKUTTA4STEPPER_HPP
