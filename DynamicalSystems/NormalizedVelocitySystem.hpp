//
// Created by Nikita Kruk on 21.11.18.
//

#ifndef SPRNEMATICALIGNMENTODEINTEGRATION_NORMALIZEDVELOCITYSYSTEM_HPP
#define SPRNEMATICALIGNMENTODEINTEGRATION_NORMALIZEDVELOCITYSYSTEM_HPP

#include "../Parallelization/Thread.hpp"
#include "../Parameters/PeriodicBoundaryConditionsConfiguration.hpp"

#include <eigen3/Eigen/Dense>

class NormalizedVelocitySystem
{
 public:

  explicit NormalizedVelocitySystem(Thread& thread, PeriodicBoundaryConditionsConfiguration &pbc_config);
  ~NormalizedVelocitySystem();

  void EvaluateRhs(std::vector<Real> &system_state,
                   const std::vector<Real> &k_prev,
                   std::vector<Real> &k_next,
                   Real k_coef,
                   Real dt);

  void EvaluateInteractionsWithAllPairs(std::vector<Real> &system_state,
                                        const std::vector<Real> &k_prev,
                                        std::vector<Real> &k_next,
                                        Real k_coef,
                                        Real dt);
  void EvaluateInteractionsWithLinkedList(std::vector<Real> &system_state,
                                          const std::vector<Real> &k_prev,
                                          std::vector<Real> &k_next,
                                          Real k_coef,
                                          Real dt);

  bool last_coefficient_;

 private:

  Thread &thread_;
  PeriodicBoundaryConditionsConfiguration &pbc_config_;

  std::vector<std::vector<Real>> alignment_force_;
  std::vector<Real> neighborhood_cardinality_;

  //Cell Subdivision Routine
  Real x_size_;
  Real y_size_;
  int num_subcells_x_;
  int num_subcells_y_;
#if defined(MPI_FAST_INTERACTION)
  std::vector<int> pre_linked_list_;
#endif
  std::vector<std::vector<int>> linked_list_;
  std::vector<std::vector<int>> neighboring_cells_;

  void AdjustNeighboringCellToPeriodicBoundaries(int &cell_x, int &cell_y);
  Eigen::Vector4d YukawaPotentialForce(Real x_i, Real y_i, Real x_j, Real y_j, Real u_x, Real u_y);
  Eigen::Vector2d YukawaForce(Real x_i, Real y_i, Real x_j, Real y_j, Real u_x, Real u_y);

};

#endif //SPRNEMATICALIGNMENTODEINTEGRATION_NORMALIZEDVELOCITYSYSTEM_HPP
