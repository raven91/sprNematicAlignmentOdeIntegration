//
// Created by Nikita Kruk on 16.11.18.
//

#ifndef SPRNEMATICALIGNMENTODEINTEGRATION_NEMATICALIGNMENTSYSTEM_HPP
#define SPRNEMATICALIGNMENTODEINTEGRATION_NEMATICALIGNMENTSYSTEM_HPP

#include "../Definitions.hpp"
#include "../Parameters/PeriodicBoundaryConditionsConfiguration.hpp"
#include "../Parallelization/Thread.hpp"

class NematicAlignmentSystem
{
 public:

  explicit NematicAlignmentSystem(Thread &thread,
                                  Real sigma,
                                  Real rho,
                                  Real D_T,
                                  Real D_R,
                                  PeriodicBoundaryConditionsConfiguration &pbc_config);
  ~NematicAlignmentSystem();

  void EvaluateRhs(std::vector<Real> &system_state,
                   const std::vector<Real> &k_prev,
                   std::vector<Real> &k_next,
                   Real k_coef,
                   Real dt);
  void EvaluateRhs(std::vector<Real> &system_state,
                   std::vector<Real> &derivative,
                   std::vector<Real> &additional_derivative,
                   Real dt);

  void EvaluateInteractionsWithAllPairs(std::vector<Real> &system_state,
                                        const std::vector<Real> &k_prev,
                                        std::vector<Real> &k_next,
                                        Real k_coef,
                                        Real dt);
  void EvaluateInteractionsWithAllPairs(std::vector<Real> &system_state,
                                        std::vector<Real> &derivative,
                                        std::vector<Real> &additional_derivative,
                                        Real dt);
  void EvaluateInteractionsWithLinkedList(std::vector<Real> &system_state,
                                          const std::vector<Real> &k_prev,
                                          std::vector<Real> &k_next,
                                          Real k_coef,
                                          Real dt);
  void EvaluateInteractionsWithLinkedList(std::vector<Real> &system_state,
                                          std::vector<Real> &derivative,
                                          std::vector<Real> &additional_derivative,
                                          Real dt);
  void EvaluateInteractionsWithVerletNeighborList(std::vector<Real> &system_state,
                                                  const std::vector<Real> &k_prev,
                                                  std::vector<Real> &k_next,
                                                  Real k_coef,
                                                  Real dt);

  Real D_T() const { return D_T_; }
  Real D_R() const { return D_R_; }

  bool last_coefficient_;

 private:

  Thread &thread_;
  PeriodicBoundaryConditionsConfiguration &pbc_config_;
  Real sigma_;
  Real rho_;
  Real rho_squared_;
  Real D_T_;
  Real D_R_;

  std::vector<Real> alignment_force_;
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

  //Verlet neighbor list
  Real r_max_;
  Real r_min_;
  std::vector<std::vector<int>> verlet_list_;
  Real accumulated_displacement_;
  bool should_update_lists_;

  void AdjustNeighboringCellToPeriodicBoundaries(int &cell_x, int &cell_y);

};

#endif //SPRNEMATICALIGNMENTODEINTEGRATION_NEMATICALIGNMENTSYSTEM_HPP
