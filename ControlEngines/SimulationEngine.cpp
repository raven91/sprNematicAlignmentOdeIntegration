//
// Created by Nikita Kruk on 16.11.18.
//

#include "SimulationEngine.hpp"

#include <cmath>
#include <fstream>
#include <cassert>
#include <algorithm> // std::copy
#include <iostream>

SimulationEngine::SimulationEngine(Thread &thread) :
    thread_(thread),
    pbc_config_(1.0, 1.0),
    system_state_(kS * kN, 0.0)
{

}

SimulationEngine::~SimulationEngine()
{
  system_state_.clear();
}

void SimulationEngine::InitializeRandomSystemState()
{
  if (thread_.IsRoot())
  {
    const Real two_pi = 2.0 * M_PI;
    std::uniform_real_distribution<Real> uni_real_dist_01(0.0, 1.0);
    std::uniform_real_distribution<Real> uni_real_dist_02pi(0.0, two_pi);
    std::uniform_real_distribution<Real> uni_real_dist_11(-1.0, 1.0);

    for (int i = 0; i < kN; ++i)
    {
      system_state_[kS * i] = uni_real_dist_01(mersenne_twister_generator);
      system_state_[kS * i + 1] = uni_real_dist_01(mersenne_twister_generator);
#if defined(CARTESIAN_REPRESENTATION)
      Real angle = uni_real_dist_02pi(mersenne_twister_generator);
      system_state_[kS * i + 2] = std::cos(angle);
      system_state_[kS * i + 3] = std::sin(angle);

//      system_state_[0] = 0.51; system_state_[1] = 0.45; system_state_[2] = 0.0; system_state_[3] = 1.0;
//      system_state_[4] = 0.5; system_state_[5] = 0.55; system_state_[6] = 0.0; system_state_[7] = -1.0;
#else
      system_state_[kS * i + 2] = uni_real_dist_02pi(mersenne_twister_generator);
#endif
    } // i

    std::cout << "system state initialization complete" << std::endl;
  }
  thread_.BroadcastVector(system_state_);
}

void SimulationEngine::RunSimulation()
{
  if (thread_.IsRoot())
  {
    std::cout << "simulation started with " << thread_.GetNumberOfMpichThreads() << " (MPICH) threads" << std::endl;
  }

  Real sigma = 0.0;
  Real rho = 0.0;
  Real D_T = 0.0;
  Real D_R = 0.0;

  Real t_0 = 0.0;
  Real t_1 = 1000.0;
  Real dt = 0.0001;

  int trial = 0;
//  for (int trial = 0; trial < 10; ++trial)
  {
    InitializeRandomSystemState();
//    StochasticRungeKuttaStepper stochastic_rk_stepper(thread_);
    RungeKutta4Stepper rk4_stepper(thread_);
//    NematicAlignmentSystem system(thread_, sigma, rho, D_T, D_R, pbc_config_);
    NormalizedVelocitySystem system(thread_, pbc_config_);
    BinaryObserver binary_observer(thread_, sigma, rho, D_T, D_R, pbc_config_, dt, trial);
//    IntegrateConst(stochastic_rk_stepper, system, t_0, t_1, dt, binary_observer);
    IntegrateConst(rk4_stepper, system, t_0, t_1, dt, binary_observer);
  } // trial

  if (thread_.IsRoot())
  {
    std::cout << "simulation ended" << std::endl;
  }
}

template<typename System>
void SimulationEngine::IntegrateConst(RungeKutta4Stepper &rk4_stepper,
                                      System &system,
                                      Real t_0,
                                      Real t_1,
                                      Real dt,
                                      BinaryObserver &binary_observer)
{
  Real t = t_0;
  binary_observer.SaveSystemState(system_state_, t);
  binary_observer.SaveSummaryStatistics(system_state_, t);
  while (t <= t_1)
  {
    t += dt;
    rk4_stepper.DoStep(system, system_state_, t, dt);
    //with the linked list technique keep all positions under periodic boundaries
    pbc_config_.ApplyPeriodicBoundaryConditions(system_state_);
    binary_observer.SaveSystemState(system_state_, t);
    binary_observer.SaveSummaryStatistics(system_state_, t);
  }
}

template<typename System>
void SimulationEngine::IntegrateConst(StochasticRungeKuttaStepper &srk_stepper,
                                      System &system,
                                      Real t_0,
                                      Real t_1,
                                      Real dt,
                                      BinaryObserver &binary_observer)
{
  Real t = t_0;
  binary_observer.SaveSystemState(system_state_, t);
  binary_observer.SaveSummaryStatistics(system_state_, t);
  while (t <= t_1)
  {
    t += dt;
    srk_stepper.DoStep(system, system_state_, t, dt);
    //with the linked list technique keep all positions under periodic boundaries
    pbc_config_.ApplyPeriodicBoundaryConditions(system_state_);
    binary_observer.SaveSystemState(system_state_, t);
    binary_observer.SaveSummaryStatistics(system_state_, t);
  }
}