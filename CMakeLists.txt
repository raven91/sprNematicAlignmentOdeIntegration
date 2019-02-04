cmake_minimum_required(VERSION 3.12)
project(sprNematicAlignmentOdeIntegration)

set(CMAKE_CXX_STANDARD 14)

add_executable(sprNematicAlignmentOdeIntegration main.cpp Observers/BinaryObserver.cpp Observers/BinaryObserver.hpp DynamicalSystems/NematicAlignmentSystem.cpp DynamicalSystems/NematicAlignmentSystem.hpp Parameters/PeriodicBoundaryConditionsConfiguration.cpp Parameters/PeriodicBoundaryConditionsConfiguration.hpp Steppers/RungeKutta4Stepper.cpp Steppers/RungeKutta4Stepper.hpp Steppers/StochasticRungeKuttaStepper.cpp Steppers/StochasticRungeKuttaStepper.hpp Definitions.hpp ControlEngines/SimulationEngine.cpp ControlEngines/SimulationEngine.hpp Parallelization/Parallelization.hpp Parallelization/Parallelization.cpp Parallelization/Thread.cpp Parallelization/Thread.hpp DynamicalSystems/NormalizedVelocitySystem.cpp DynamicalSystems/NormalizedVelocitySystem.hpp)

include_directories(/Users/nikita/Documents/boost/bin/include/ /Users/nikita/Documents/mpich-install/include/)
link_directories(/Users/nikita/Documents/boost/bin/lib/ /Users/nikita/Documents/mpich-install/lib/)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()
CHECK_CXX_COMPILER_FLAG("-O3" COMPILER_SUPPORTS_O3)
if(COMPILER_SUPPORTS_O3)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif()
CHECK_CXX_COMPILER_FLAG("-flto" COMPILER_SUPPORTS_FLTO)
if(COMPILER_SUPPORTS_FLTO)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
endif()