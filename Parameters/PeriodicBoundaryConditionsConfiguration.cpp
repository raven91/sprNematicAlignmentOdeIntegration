//
// Created by Nikita Kruk on 16.11.18.
//

#include "PeriodicBoundaryConditionsConfiguration.hpp"

PeriodicBoundaryConditionsConfiguration::PeriodicBoundaryConditionsConfiguration(Real x_size, Real y_size) :
    x_size_(x_size),
    y_size_(y_size),
    x_rsize_(1.0 / x_size),
    y_rsize_(1.0 / y_size)
{

}

PeriodicBoundaryConditionsConfiguration::~PeriodicBoundaryConditionsConfiguration()
{

}

// if all interaction sites are in the same simulation box
void PeriodicBoundaryConditionsConfiguration::ClassAEffectiveParticleDistance(Real x_i,
                                                                              Real y_i,
                                                                              Real x_j,
                                                                              Real y_j,
                                                                              Real &dx,
                                                                              Real &dy)
{
  dx = x_j - x_i;
  dx -= static_cast<int>(dx * 2.0 * x_rsize_) * x_size_;

  dy = y_j - y_i;
  dy -= static_cast<int>(dy * 2.0 * y_rsize_) * y_size_;
}

// if the centers or primary sites of all molecules are in the same box, but the other sites can be away from the centers
void PeriodicBoundaryConditionsConfiguration::ClassBEffectiveParticleDistanceSigned(Real x_i,
                                                                                    Real y_i,
                                                                                    Real x_j,
                                                                                    Real y_j,
                                                                                    Real &dx,
                                                                                    Real &dy)
{
  dx = x_j - x_i;
  if (dx > 0.5 * x_size_)
  {
    dx -= x_size_;
  } else if (dx < -0.5 * x_size_)
  {
    dx += x_size_;
  }

  dy = y_j - y_i;
  if (dy > 0.5 * y_size_)
  {
    dy -= y_size_;
  } else if (dy < -0.5 * y_size_)
  {
    dy += y_size_;
  }
}

// if the centers or primary sites of all molecules are in the same box, but the other sites can be away from the centers
void PeriodicBoundaryConditionsConfiguration::ClassBEffectiveParticleDistanceUnsigned(Real x_i,
                                                                                      Real y_i,
                                                                                      Real x_j,
                                                                                      Real y_j,
                                                                                      Real &dx,
                                                                                      Real &dy)
{
  dx = std::fabs(x_j - x_i);
  if (dx > 0.5 * x_size_)
  {
    dx -= x_size_;
  }

  dy = std::fabs(y_j - y_i);
  if (dy > 0.5 * y_size_)
  {
    dy -= y_size_;
  }
}

// calculate remaindars at any distance
//if the sign of the distance is relevant
void PeriodicBoundaryConditionsConfiguration::ClassCEffectiveParticleDistanceSigned(Real x_i,
                                                                                    Real y_i,
                                                                                    Real x_j,
                                                                                    Real y_j,
                                                                                    Real &dx,
                                                                                    Real &dy)
{
  dx = x_j - x_i;
  dx -= x_size_ * std::nearbyint(dx * x_rsize_);

  dy = y_j - y_i;
  dy -= y_size_ * std::nearbyint(dy * y_rsize_);
}

// calculate remaindars at any distance
//if the sign of the distance is not relevant
//#pragma acc routine seq
void PeriodicBoundaryConditionsConfiguration::ClassCEffectiveParticleDistanceUnsigned(Real x_i,
                                                                                      Real y_i,
                                                                                      Real x_j,
                                                                                      Real y_j,
                                                                                      Real &dx,
                                                                                      Real &dy)
{
  dx = std::fabs(x_j - x_i);
  dx -= static_cast<int>(dx * x_rsize_ + 0.5) * x_size_;

  dy = std::fabs(y_j - y_i);
  dy -= static_cast<int>(dy * y_rsize_ + 0.5) * y_size_;
}

void PeriodicBoundaryConditionsConfiguration::ApplyPeriodicBoundaryConditions(Real x, Real y, Real &x_pbc, Real &y_pbc)
{
  x_pbc = x - std::floor(x * x_rsize_) * x_size_;
  y_pbc = y - std::floor(y * y_rsize_) * y_size_;
}

void PeriodicBoundaryConditionsConfiguration::ApplyPeriodicBoundaryConditions(std::vector<Real> &system_state)
{
#if defined(MPI_FAST_INTERACTION)
  for (int i = mpi_thread_rank * mpi_num_particles_per_thread; i < (mpi_thread_rank + 1) * mpi_num_particles_per_thread; ++i)
#else
  for (int i = 0; i < kN; ++i)
#endif
  {
    system_state[kS * i] -= std::floor(system_state[kS * i] * x_rsize_) * x_size_;
    system_state[kS * i + 1] -= std::floor(system_state[kS * i + 1] * y_rsize_) * y_size_;
  }
}