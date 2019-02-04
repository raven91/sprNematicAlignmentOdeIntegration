#include "ControlEngines/SimulationEngine.hpp"
#include "Parallelization/Parallelization.hpp"
#include "Parallelization/Thread.hpp"

// initialization of the random number generator
std::mt19937 mersenne_twister_generator(std::random_device{}());

int main(int argc, char **argv)
{
  LaunchParallelSession(argc, argv);

  Thread thread(argc, argv);
  SimulationEngine engine(thread);
  engine.RunSimulation();

  FinalizeParallelSession();

  return 0;
}