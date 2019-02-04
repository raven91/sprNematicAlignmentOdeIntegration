//
// Created by Nikita Kruk on 16.11.18.
//

#include "RungeKutta4Stepper.hpp"

#include <algorithm> // std::copy
#include <iterator>  // std::advance

RungeKutta4Stepper::RungeKutta4Stepper(Thread &thread) :
    thread_(thread),
    k_1_(kS * kN, 0.0),
    k_2_(kS * kN, 0.0),
    k_3_(kS * kN, 0.0),
    k_4_(kS * kN, 0.0)
{

}

RungeKutta4Stepper::~RungeKutta4Stepper()
{
  k_1_.clear();
  k_2_.clear();
  k_3_.clear();
  k_4_.clear();
}

void RungeKutta4Stepper::DoStep(NematicAlignmentSystem &system, std::vector<Real> &system_state, Real t, Real dt)
{
  // k_1
//	chimera_system.first_coefficient_ = true;
  // No need to reset k_1_ to 0.0, since k_coef is 0.0
  system.EvaluateRhs(system_state, k_1_, k_1_, 0.0, dt);
  thread_.SynchronizeVectorThroughBuffer(k_1_);
//	chimera_system.first_coefficient_ = false;

  // k_2
  system.EvaluateRhs(system_state, k_1_, k_2_, 0.5, dt);
  thread_.SynchronizeVectorThroughBuffer(k_2_);

  // k_3
  system.EvaluateRhs(system_state, k_2_, k_3_, 0.5, dt);
  thread_.SynchronizeVectorThroughBuffer(k_3_);

  // k_4
  system.last_coefficient_ = true;
  system.EvaluateRhs(system_state, k_3_, k_4_, 1.0, dt);
  system.last_coefficient_ = false;

  for (int i : thread_.GetLoopIndices())
  {
    int ii = kS * i;
    system_state[ii] += (k_1_[ii] + 2.0 * k_2_[ii] + 2.0 * k_3_[ii] + k_4_[ii]) * dt / 6.0;
    system_state[ii + 1] += (k_1_[ii + 1] + 2.0 * k_2_[ii + 1] + 2.0 * k_3_[ii + 1] + k_4_[ii + 1]) * dt / 6.0;
    system_state[ii + 2] += (k_1_[ii + 2] + 2.0 * k_2_[ii + 2] + 2.0 * k_3_[ii + 2] + k_4_[ii + 2]) * dt / 6.0;
#if defined(CARTESIAN_REPRESENTATION)
    system_state[ii + 3] += (k_1_[ii + 3] + 2.0 * k_2_[ii + 3] + 2.0 * k_3_[ii + 3] + k_4_[ii + 3]) * dt / 6.0;
#endif
  } // i
  thread_.SynchronizeVectorThroughBuffer(system_state);
}

void RungeKutta4Stepper::DoStep(NormalizedVelocitySystem &system, std::vector<Real> &system_state, Real t, Real dt)
{
  // k_1
//	chimera_system.first_coefficient_ = true;
  // No need to reset k_1_ to 0.0, since k_coef is 0.0
  system.EvaluateRhs(system_state, k_1_, k_1_, 0.0, dt);
  thread_.SynchronizeVectorThroughBuffer(k_1_);
//	chimera_system.first_coefficient_ = false;

  // k_2
  system.EvaluateRhs(system_state, k_1_, k_2_, 0.5, dt);
  thread_.SynchronizeVectorThroughBuffer(k_2_);

  // k_3
  system.EvaluateRhs(system_state, k_2_, k_3_, 0.5, dt);
  thread_.SynchronizeVectorThroughBuffer(k_3_);

  // k_4
  system.last_coefficient_ = true;
  system.EvaluateRhs(system_state, k_3_, k_4_, 1.0, dt);
  system.last_coefficient_ = false;

  for (int i : thread_.GetLoopIndices())
  {
    int ii = kS * i;
    system_state[ii] += (k_1_[ii] + 2.0 * k_2_[ii] + 2.0 * k_3_[ii] + k_4_[ii]) * dt / 6.0;
    system_state[ii + 1] += (k_1_[ii + 1] + 2.0 * k_2_[ii + 1] + 2.0 * k_3_[ii + 1] + k_4_[ii + 1]) * dt / 6.0;
    system_state[ii + 2] += (k_1_[ii + 2] + 2.0 * k_2_[ii + 2] + 2.0 * k_3_[ii + 2] + k_4_[ii + 2]) * dt / 6.0;
#if defined(CARTESIAN_REPRESENTATION)
    system_state[ii + 3] += (k_1_[ii + 3] + 2.0 * k_2_[ii + 3] + 2.0 * k_3_[ii + 3] + k_4_[ii + 3]) * dt / 6.0;
#endif
  } // i
  thread_.SynchronizeVectorThroughBuffer(system_state);
}