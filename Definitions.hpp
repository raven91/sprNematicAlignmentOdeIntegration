//
// Created by Nikita Kruk on 16.11.18.
//

#ifndef SPRNEMATICALIGNMENTODEINTEGRATION_DEFINITIONS_HPP
#define SPRNEMATICALIGNMENTODEINTEGRATION_DEFINITIONS_HPP

//#define LICHTENBERG
//#define BCS_CLUSTER

#define CARTESIAN_REPRESENTATION

#include <random>

typedef double Real;

const int kN = 100;
#if defined(CARTESIAN_REPRESENTATION)
const int kS = 4;
#else
const int kS = 3;
#endif

extern std::mt19937 mersenne_twister_generator;

#endif //SPRNEMATICALIGNMENTODEINTEGRATION_DEFINITIONS_HPP
