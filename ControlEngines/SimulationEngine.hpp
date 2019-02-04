//
// Created by Nikita Kruk on 16.11.18.
//

#ifndef SPRNEMATICALIGNMENTODEINTEGRATION_SIMULATIONENGINE_HPP
#define SPRNEMATICALIGNMENTODEINTEGRATION_SIMULATIONENGINE_HPP

#include "../Definitions.hpp"
#include "../Parameters/PeriodicBoundaryConditionsConfiguration.hpp"
#include "../Steppers/RungeKutta4Stepper.hpp"
#include "../Steppers/StochasticRungeKuttaStepper.hpp"
#include "../DynamicalSystems/NematicAlignmentSystem.hpp"
#include "../DynamicalSystems/NormalizedVelocitySystem.hpp"
#include "../Observers/BinaryObserver.hpp"
#include "../Parallelization/Thread.hpp"

#include <string>

class SimulationEngine
{
 public:

  explicit SimulationEngine(Thread &thread);
  ~SimulationEngine();

  void RunSimulation();

 private:

  Thread &thread_;
  PeriodicBoundaryConditionsConfiguration pbc_config_;
  std::vector<Real> system_state_;

  void InitializeRandomSystemState();
  void InitializeSystemStateFromFile(const std::string &file_name);

  template<typename System>
  void IntegrateConst(RungeKutta4Stepper &rk4_stepper,
                      System &system,
                      Real t_0,
                      Real t_1,
                      Real dt,
                      BinaryObserver &binary_observer);
  template<typename System>
  void IntegrateConst(StochasticRungeKuttaStepper &srk_stepper,
                      System &system,
                      Real t_0,
                      Real t_1,
                      Real dt,
                      BinaryObserver &binary_observer);

};

#endif //SPRNEMATICALIGNMENTODEINTEGRATION_SIMULATIONENGINE_HPP
