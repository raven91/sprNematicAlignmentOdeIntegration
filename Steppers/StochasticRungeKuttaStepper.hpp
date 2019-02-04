//
// Created by Nikita Kruk on 16.11.18.
//

#ifndef SPRNEMATICALIGNMENTODEINTEGRATION_STOCHASTICRUNGEKUTTASTEPPER_HPP
#define SPRNEMATICALIGNMENTODEINTEGRATION_STOCHASTICRUNGEKUTTASTEPPER_HPP

#include "../Definitions.hpp"
#include "../DynamicalSystems/NematicAlignmentSystem.hpp"
#include "../Parallelization/Thread.hpp"

#include <cmath>
#include <vector>
#include <random>

/**
 * Numerical Solution of Stochastic Differential Equations with Jumps in Finance
 * Strong Order 1.5 Taylor Scheme
 */
class StochasticRungeKuttaStepper
{
 public:

  explicit StochasticRungeKuttaStepper(Thread &thread) : thread_(thread)
  {};
  ~StochasticRungeKuttaStepper()
  {};

  void DoStep(NematicAlignmentSystem &system, std::vector<Real> &system_state, Real t, Real dt)
  {
    std::vector<Real> derivative(kS * kN, 0.0), additional_derivative(kN, 0.0);
    system.EvaluateRhs(system_state, derivative, additional_derivative, dt);

    std::normal_distribution<Real> norm_dist(0.0, 1.0);
    Real U_1 = 0.0, U_2 = 0.0, delta_W = 0.0, delta_Z = 0.0;

    for (int i : thread_.GetLoopIndices())
    {
      U_1 = norm_dist(mersenne_twister_generator);
      U_2 = norm_dist(mersenne_twister_generator);
      delta_W = U_1 * std::sqrt(dt);
      delta_Z = 0.5 * dt * std::sqrt(dt) * (U_1 + U_2 / std::sqrt(3.0));

      int ii = kS * i;
      system_state[ii] += derivative[ii] * dt - std::sqrt(2.0 * system.D_R()) * derivative[ii + 1] * delta_Z
          - 0.5 * (derivative[ii + 2] * derivative[ii + 1] + system.D_R() * derivative[ii]) * dt * dt;
      system_state[ii + 1] += derivative[ii + 1] * dt + std::sqrt(2.0 * system.D_R()) * derivative[ii] * delta_Z
          + 0.5 * (derivative[ii + 2] * derivative[ii] - system.D_R() * derivative[ii + 1]) * dt * dt;
      system_state[ii + 2] += derivative[ii + 2] * dt + std::sqrt(2.0 * system.D_R()) * delta_W
          - std::sqrt(2.0 * system.D_R()) * additional_derivative[i] * delta_Z
          - 0.5 * (derivative[ii + 2] * additional_derivative[i] + system.D_R() * derivative[ii + 2]) * dt * dt;
    } // i
    thread_.SynchronizeVectorThroughBuffer(system_state);
  }

 private:

  Thread &thread_;

};

#endif //SPRNEMATICALIGNMENTODEINTEGRATION_STOCHASTICRUNGEKUTTASTEPPER_HPP
