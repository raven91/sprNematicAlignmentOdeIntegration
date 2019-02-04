//
// Created by Nikita Kruk on 16.11.18.
//

#ifndef SPRNEMATICALIGNMENTODEINTEGRATION_PERIODICBOUNDARYCONDITIONSCONFIGURATION_HPP
#define SPRNEMATICALIGNMENTODEINTEGRATION_PERIODICBOUNDARYCONDITIONSCONFIGURATION_HPP

#include "../Definitions.hpp"

class PeriodicBoundaryConditionsConfiguration
{
 public:

  PeriodicBoundaryConditionsConfiguration(Real x_size, Real y_size);
  ~PeriodicBoundaryConditionsConfiguration();

  void ClassAEffectiveParticleDistance(Real x_i, Real y_i, Real x_j, Real y_j, Real &dx, Real &dy);

  void ClassBEffectiveParticleDistanceSigned(Real x_i, Real y_i, Real x_j, Real y_j, Real &dx, Real &dy);
  void ClassBEffectiveParticleDistanceUnsigned(Real x_i, Real y_i, Real x_j, Real y_j, Real &dx, Real &dy);

  void ClassCEffectiveParticleDistanceSigned(Real x_i, Real y_i, Real x_j, Real y_j, Real &dx, Real &dy);
  void ClassCEffectiveParticleDistanceUnsigned(Real x_i, Real y_i, Real x_j, Real y_j, Real &dx, Real &dy);

  void ApplyPeriodicBoundaryConditions(Real x, Real y, Real &x_pbc, Real &y_pbc);
  void ApplyPeriodicBoundaryConditions(std::vector<Real> &system_state);

 private:

  Real x_size_;
  Real y_size_;
  Real x_rsize_;
  Real y_rsize_;

};

#endif //SPRNEMATICALIGNMENTODEINTEGRATION_PERIODICBOUNDARYCONDITIONSCONFIGURATION_HPP
