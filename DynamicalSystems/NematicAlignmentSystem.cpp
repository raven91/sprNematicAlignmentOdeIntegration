//
// Created by Nikita Kruk on 16.11.18.
//

#include "NematicAlignmentSystem.hpp"

#include <cmath>
#include <algorithm> // std::fill

NematicAlignmentSystem::NematicAlignmentSystem(Thread &thread,
                                               Real sigma,
                                               Real rho,
                                               Real D_T,
                                               Real D_R,
                                               PeriodicBoundaryConditionsConfiguration &pbc_config) :
    thread_(thread),
    pbc_config_(pbc_config),
    sigma_(sigma),
    rho_(rho),
    rho_squared_(rho * rho),
    D_T_(D_T),
    D_R_(D_R),
    alignment_force_(kN, 0.0),
    neighborhood_cardinality_(kN, 0.0),
    x_size_(1.0),
    y_size_(1.0),
    num_subcells_x_(1.0 / rho),
    num_subcells_y_(1.0 / rho)
{
#if defined(MPI_FAST_INTERACTION)
  pre_linked_list_ = std::vector<int>(kN, 0);
#endif

  linked_list_ = std::vector<std::vector<int>>(num_subcells_x_ * num_subcells_y_, std::vector<int>());

  //all neighbors
  neighboring_cells_ =
      {
          {-1, -1}, {0, -1}, {1, -1},
          {-1, 0}, {0, 0}, {1, 0},
          {-1, 1}, {0, 1}, {1, 1}
      };
  //half of neighbors
//	neighboring_cells_ =
//	{
//			{0, 0}, {1, 0},
//			{-1, 1}, {0, 1}, {1, 1}
//	};

//Verlet neighbor list
  last_coefficient_ = true;
  r_min_ = rho;
  r_max_ = r_min_ + 2.0 * 1 * 0.01;//2 * k * dt
  should_update_lists_ = true;
  accumulated_displacement_ = 0.0;
#if defined(MPI_FAST_INTERACTION)
  verlet_list_ = std::vector<std::vector<int>>(mpi_num_particles_per_thread, std::vector<int>());
#else
  verlet_list_ = std::vector<std::vector<int>>(kN, std::vector<int>());
#endif
}

NematicAlignmentSystem::~NematicAlignmentSystem()
{
  alignment_force_.clear();
  neighborhood_cardinality_.clear();
  neighboring_cells_.clear();
  for (int i = 0; i < linked_list_.size(); ++i)
  {
    linked_list_[i].clear();
  }
  linked_list_.clear();
  for (int i = 0; i < verlet_list_.size(); ++i)
  {
    verlet_list_[i].clear();
  }
  verlet_list_.clear();
}

void NematicAlignmentSystem::EvaluateRhs(std::vector<Real> &system_state,
                                         const std::vector<Real> &k_prev,
                                         std::vector<Real> &k_next,
                                         Real k_coef,
                                         Real dt)
{
//#ifdef LINKED_LIST
//	if (rho_ <= 0.01)
//	{
//		EvaluateInteractionsWithVerletNeighborList(system_state, k_prev, k_next, k_coef, dt);
//	}
//	else
  if (rho_ <= 0.25)
  {
    EvaluateInteractionsWithLinkedList(system_state, k_prev, k_next, k_coef, dt);
  }
//#else
  else
  {
    EvaluateInteractionsWithAllPairs(system_state, k_prev, k_next, k_coef, dt);
  }
//#endif
}

void NematicAlignmentSystem::EvaluateRhs(std::vector<Real> &system_state,
                                         std::vector<Real> &derivative,
                                         std::vector<Real> &additional_derivative,
                                         Real dt)
{
  if (rho_ <= 0.25)
  {
    EvaluateInteractionsWithLinkedList(system_state, derivative, additional_derivative, dt);
  } else
  {
    EvaluateInteractionsWithAllPairs(system_state, derivative, additional_derivative, dt);
  }
}

//Calculate the deterministic part of ODE/SDE
void NematicAlignmentSystem::EvaluateInteractionsWithAllPairs(std::vector<Real> &system_state,
                                                              const std::vector<Real> &k_prev,
                                                              std::vector<Real> &k_next,
                                                              Real k_coef,
                                                              Real dt)
{
  std::fill(alignment_force_.begin(), alignment_force_.end(), 0.0);
  std::fill(neighborhood_cardinality_.begin(), neighborhood_cardinality_.end(), 0.0);

  Real x_i = 0.0, y_i = 0.0, phi_i = 0.0;
  Real x_j = 0.0, y_j = 0.0, phi_j = 0.0;
  Real dx = 0.0, dy = 0.0, dr_squared = 0.0;

  for (int i : thread_.GetLoopIndices())
  {
    int ii = kS * i;

    x_i = system_state[ii] + k_coef * dt * k_prev[ii];
    y_i = system_state[ii + 1] + k_coef * dt * k_prev[ii + 1];
    phi_i = system_state[ii + 2] + k_coef * dt * k_prev[ii + 2];

    //j ~ i
    for (int j = 0; j < kN; ++j)
    {
      int jj = kS * j;

      if (j != i)
      {
        x_j = system_state[jj] + k_coef * dt * k_prev[jj];
        y_j = system_state[jj + 1] + k_coef * dt * k_prev[jj + 1];
        phi_j = system_state[jj + 2] + k_coef * dt * k_prev[jj + 2];

        pbc_config_.ClassAEffectiveParticleDistance(x_i, y_i, x_j, y_j, dx, dy);
        dr_squared = dx * dx + dy * dy;

        if (dr_squared <= rho_squared_)
        {
          alignment_force_[i] += std::sin(2.0 * (phi_j - phi_i));
          ++neighborhood_cardinality_[i];
        }
      }
    } // j

    k_next[ii] = std::cos(phi_i);
    k_next[ii + 1] = std::sin(phi_i);
    if (neighborhood_cardinality_[i] != 0.0)
    {
      k_next[ii + 2] = sigma_ * alignment_force_[i] / neighborhood_cardinality_[i];
    }
  } // i
}

void NematicAlignmentSystem::EvaluateInteractionsWithAllPairs(std::vector<Real> &system_state,
                                                              std::vector<Real> &derivative,
                                                              std::vector<Real> &additional_derivative,
                                                              Real dt)
{
  std::fill(alignment_force_.begin(), alignment_force_.end(), 0.0);
  std::fill(neighborhood_cardinality_.begin(), neighborhood_cardinality_.end(), 0.0);

  Real x_i = 0.0, y_i = 0.0, phi_i = 0.0;
  Real x_j = 0.0, y_j = 0.0, phi_j = 0.0;
  Real dx = 0.0, dy = 0.0, dr_squared = 0.0;

  for (int i : thread_.GetLoopIndices())
  {
    int ii = kS * i;

    x_i = system_state[ii];
    y_i = system_state[ii + 1];
    phi_i = system_state[ii + 2];

    //j ~ i
    for (int j = 0; j < kN; ++j)
    {
      int jj = kS * j;

      if (j != i)
      {
        x_j = system_state[jj];
        y_j = system_state[jj + 1];
        phi_j = system_state[jj + 2];

        pbc_config_.ClassAEffectiveParticleDistance(x_i, y_i, x_j, y_j, dx, dy);
        dr_squared = dx * dx + dy * dy;

        if (dr_squared <= rho_squared_)
        {
          alignment_force_[i] += std::sin(2.0 * (phi_j - phi_i));
          ++neighborhood_cardinality_[i];
          additional_derivative[i] += 2.0 * std::cos(2.0 * (phi_j - phi_i));
        }
      }
    } // j

    derivative[ii] = std::cos(phi_i);
    derivative[ii + 1] = std::sin(phi_i);
    if (neighborhood_cardinality_[i] != 0.0)
    {
      derivative[ii + 2] = sigma_ * alignment_force_[i] / neighborhood_cardinality_[i];
      additional_derivative[i] *= sigma_ / neighborhood_cardinality_[i];
    }
  }
}

//Calculate the deterministic part of ODE/SDE
void NematicAlignmentSystem::EvaluateInteractionsWithLinkedList(std::vector<Real> &system_state,
                                                                const std::vector<Real> &k_prev,
                                                                std::vector<Real> &k_next,
                                                                Real k_coef,
                                                                Real dt)
{
  std::fill(alignment_force_.begin(), alignment_force_.end(), 0.0);
  std::fill(neighborhood_cardinality_.begin(), neighborhood_cardinality_.end(), 0.0);

  Real x_i = 0.0, y_i = 0.0, phi_i = 0.0;
  Real x_j = 0.0, y_j = 0.0, phi_j = 0.0;
  Real dx = 0.0, dy = 0.0, dr_squared = 0.0;

  for (int i = 0; i < linked_list_.size(); ++i)
  {
    linked_list_[i].clear();
  }

  // construct a linked list
  int i_cell = 0, i_cell_x = 0, i_cell_y = 0;
  int ii = 0;
//	if (first_coefficient_)
//	{
#if defined(MPI_FAST_INTERACTION)
  for (int i : thread_.GetLoopIndices())
    {
        ii = kS * i;

        x_i   = system_state[ii    ] + k_coef * dt * k_prev[ii    ];
        y_i   = system_state[ii + 1] + k_coef * dt * k_prev[ii + 1];
//		phi_i = system_state[ii + 2] + k_coef * dt * k_prev[ii + 2];

        if (x_i < 0.0 || x_i >= x_size_ || y_i < 0.0 || y_i >= y_size_)
        {
            pbc_config_.ApplyPeriodicBoundaryConditions(x_i, y_i, x_i, y_i);
        }

        i_cell_x = int(x_i / x_size_ * num_subcells_x_);
        i_cell_y = int(y_i / y_size_ * num_subcells_y_);
        i_cell = i_cell_y * num_subcells_x_ + i_cell_x;

        pre_linked_list_[i] = i_cell;
    }

    if (mpi_thread_rank != mpi_root_rank)
    {
        MPI_Gather(&pre_linked_list_[mpi_thread_rank * mpi_num_particles_per_thread], mpi_num_particles_per_thread, MPI_INT, &pre_linked_list_[0], mpi_num_particles_per_thread, MPI_INT, mpi_root_rank, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Gather(MPI_IN_PLACE, mpi_num_particles_per_thread, MPI_INT, &pre_linked_list_[mpi_root_rank * mpi_num_particles_per_thread], mpi_num_particles_per_thread, MPI_INT, mpi_root_rank, MPI_COMM_WORLD);
    }
    MPI_Bcast(&pre_linked_list_[0], (int)pre_linked_list_.size(), MPI_INT, mpi_root_rank, MPI_COMM_WORLD);

    for (int i = 0; i < kN; ++i)
    {
        linked_list_[pre_linked_list_[i]].push_back(i);
    }
#else
  for (int i : thread_.GetLoopIndices())
  {
    ii = kS * i;

    x_i = system_state[ii] + k_coef * dt * k_prev[ii];
    y_i = system_state[ii + 1] + k_coef * dt * k_prev[ii + 1];
//		phi_i = system_state[ii + 2] + k_coef * dt * k_prev[ii + 2];

    if (x_i < 0.0 || x_i >= x_size_ || y_i < 0.0 || y_i >= y_size_)
    {
      pbc_config_.ApplyPeriodicBoundaryConditions(x_i, y_i, x_i, y_i);
    }

    i_cell_x = int(x_i / x_size_ * num_subcells_x_);
    i_cell_y = int(y_i / y_size_ * num_subcells_y_);
    i_cell = i_cell_y * num_subcells_x_ + i_cell_x;

    linked_list_[i_cell].push_back(i);
  }
#endif
//	}

  //loop using the linked list
  int j = 0, jj = 0;
  int j_cell = 0, j_cell_x = 0, j_cell_y = 0;
  for (int i : thread_.GetLoopIndices())
  {
    ii = kS * i;

    x_i = system_state[ii] + k_coef * dt * k_prev[ii];
    y_i = system_state[ii + 1] + k_coef * dt * k_prev[ii + 1];
    phi_i = system_state[ii + 2] + k_coef * dt * k_prev[ii + 2];

    if (x_i < 0.0 || x_i >= x_size_ || y_i < 0.0 || y_i >= y_size_)
    {
      pbc_config_.ApplyPeriodicBoundaryConditions(x_i, y_i, x_i, y_i);
    }

    i_cell_x = int(x_i / x_size_ * num_subcells_x_);
    i_cell_y = int(y_i / y_size_ * num_subcells_y_);
    i_cell = i_cell_y * num_subcells_x_ + i_cell_x;

    for (int neighboring_cell = 0; neighboring_cell < neighboring_cells_.size(); ++neighboring_cell)
    {
      j_cell_x = i_cell_x + neighboring_cells_[neighboring_cell][0];
      j_cell_y = i_cell_y + neighboring_cells_[neighboring_cell][1];
      AdjustNeighboringCellToPeriodicBoundaries(j_cell_x, j_cell_y);
      j_cell = j_cell_y * num_subcells_x_ + j_cell_x;

      for (std::vector<int>::iterator neighbr_iter = linked_list_[j_cell].begin();
           neighbr_iter != linked_list_[j_cell].end(); ++neighbr_iter)
      {
        j = *neighbr_iter;
        jj = kS * j;

        if (j != i)
        {
          x_j = system_state[jj] + k_coef * dt * k_prev[jj];
          y_j = system_state[jj + 1] + k_coef * dt * k_prev[jj + 1];
          phi_j = system_state[jj + 2] + k_coef * dt * k_prev[jj + 2];

          if (x_j < 0.0 || x_j >= x_size_ || y_j < 0.0 || y_j >= y_size_)
          {
            pbc_config_.ApplyPeriodicBoundaryConditions(x_j, y_j, x_j, y_j);
          }

          pbc_config_.ClassAEffectiveParticleDistance(x_i, y_i, x_j, y_j, dx, dy);
          dr_squared = dx * dx + dy * dy;

          if (dr_squared <= rho_squared_)
          {
            alignment_force_[i] += std::sin(2.0 * (phi_j - phi_i));
            ++neighborhood_cardinality_[i];
          }
        }
      } // neighbr_iter
    } // neighboring_cell

    k_next[ii] = std::cos(phi_i);
    k_next[ii + 1] = std::sin(phi_i);
    if (neighborhood_cardinality_[i] != 0)
    {
      k_next[ii + 2] = sigma_ * alignment_force_[i] / neighborhood_cardinality_[i];
    }
  } // i
}

void NematicAlignmentSystem::EvaluateInteractionsWithLinkedList(std::vector<Real> &system_state,
                                                                std::vector<Real> &derivative,
                                                                std::vector<Real> &additional_derivative,
                                                                Real dt)
{
  std::fill(alignment_force_.begin(), alignment_force_.end(), 0.0);
  std::fill(neighborhood_cardinality_.begin(), neighborhood_cardinality_.end(), 0.0);

  Real x_i = 0.0, y_i = 0.0, phi_i = 0.0;
  Real x_j = 0.0, y_j = 0.0, phi_j = 0.0;
  Real dx = 0.0, dy = 0.0, dr_squared = 0.0;

  for (int i = 0; i < linked_list_.size(); ++i)
  {
    linked_list_[i].clear();
  }

  // construct a linked list
  int i_cell = 0, i_cell_x = 0, i_cell_y = 0;
  int ii = 0;
  //	if (first_coefficient_)
  //	{
#if defined(MPI_FAST_INTERACTION)
  for (int i : thread_.GetLoopIndices())
    {
        ii = kS * i;

        x_i   = system_state[ii    ];
        y_i   = system_state[ii + 1];
        //		phi_i = system_state[ii + 2];

        if (x_i < 0.0 || x_i >= x_size_ || y_i < 0.0 || y_i >= y_size_)
        {
            pbc_config_.ApplyPeriodicBoundaryConditions(x_i, y_i, x_i, y_i);
        }

        i_cell_x = int(x_i / x_size_ * num_subcells_x_);
        i_cell_y = int(y_i / y_size_ * num_subcells_y_);
        i_cell = i_cell_y * num_subcells_x_ + i_cell_x;

        pre_linked_list_[i] = i_cell;
    }

    if (mpi_thread_rank != mpi_root_rank)
    {
        MPI_Gather(&pre_linked_list_[mpi_thread_rank * mpi_num_particles_per_thread], mpi_num_particles_per_thread, MPI_INT, &pre_linked_list_[0], mpi_num_particles_per_thread, MPI_INT, mpi_root_rank, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Gather(MPI_IN_PLACE, mpi_num_particles_per_thread, MPI_INT, &pre_linked_list_[mpi_root_rank * mpi_num_particles_per_thread], mpi_num_particles_per_thread, MPI_INT, mpi_root_rank, MPI_COMM_WORLD);
    }
    MPI_Bcast(&pre_linked_list_[0], (int)pre_linked_list_.size(), MPI_INT, mpi_root_rank, MPI_COMM_WORLD);

    for (int i = 0; i < kN; ++i)
    {
        linked_list_[pre_linked_list_[i]].push_back(i);
    }
#else
  for (int i : thread_.GetLoopIndices())
  {
    ii = kS * i;

    x_i = system_state[ii];
    y_i = system_state[ii + 1];
    //		phi_i = system_state[ii + 2];

    if (x_i < 0.0 || x_i >= x_size_ || y_i < 0.0 || y_i >= y_size_)
    {
      pbc_config_.ApplyPeriodicBoundaryConditions(x_i, y_i, x_i, y_i);
    }

    i_cell_x = int(x_i / x_size_ * num_subcells_x_);
    i_cell_y = int(y_i / y_size_ * num_subcells_y_);
    i_cell = i_cell_y * num_subcells_x_ + i_cell_x;

    linked_list_[i_cell].push_back(i);
  }
#endif
  //	}

  //loop using the linked list
  int j = 0, jj = 0;
  int j_cell = 0, j_cell_x = 0, j_cell_y = 0;
  for (int i : thread_.GetLoopIndices())
  {
    ii = kS * i;

    x_i = system_state[ii];
    y_i = system_state[ii + 1];
    phi_i = system_state[ii + 2];

    if (x_i < 0.0 || x_i >= x_size_ || y_i < 0.0 || y_i >= y_size_)
    {
      pbc_config_.ApplyPeriodicBoundaryConditions(x_i, y_i, x_i, y_i);
    }

    i_cell_x = int(x_i / x_size_ * num_subcells_x_);
    i_cell_y = int(y_i / y_size_ * num_subcells_y_);
    i_cell = i_cell_y * num_subcells_x_ + i_cell_x;

    for (int neighboring_cell = 0; neighboring_cell < neighboring_cells_.size(); ++neighboring_cell)
    {
      j_cell_x = i_cell_x + neighboring_cells_[neighboring_cell][0];
      j_cell_y = i_cell_y + neighboring_cells_[neighboring_cell][1];
      AdjustNeighboringCellToPeriodicBoundaries(j_cell_x, j_cell_y);
      j_cell = j_cell_y * num_subcells_x_ + j_cell_x;

      for (std::vector<int>::iterator neighbr_iter = linked_list_[j_cell].begin();
           neighbr_iter != linked_list_[j_cell].end(); ++neighbr_iter)
      {
        j = *neighbr_iter;
        jj = kS * j;

        if (j != i)
        {
          x_j = system_state[jj];
          y_j = system_state[jj + 1];
          phi_j = system_state[jj + 2];

          if (x_j < 0.0 || x_j >= x_size_ || y_j < 0.0 || y_j >= y_size_)
          {
            pbc_config_.ApplyPeriodicBoundaryConditions(x_j, y_j, x_j, y_j);
          }

          pbc_config_.ClassAEffectiveParticleDistance(x_i, y_i, x_j, y_j, dx, dy);
          dr_squared = dx * dx + dy * dy;

          if (dr_squared <= rho_squared_)
          {
            alignment_force_[i] += std::sin(2.0 * (phi_j - phi_i));
            ++neighborhood_cardinality_[i];
            additional_derivative[i] += 2.0 * std::cos(2.0 * (phi_j - phi_i));
          }
        }
      } // neighbr_iter
    } // neighboring_cell

    derivative[ii] = 0.5 * std::cos(phi_i);
    derivative[ii + 1] = 0.5 * std::sin(phi_i);
    if (neighborhood_cardinality_[i] != 0)
    {
      derivative[ii + 2] = sigma_ * alignment_force_[i] / neighborhood_cardinality_[i];
      additional_derivative[i] *= sigma_ / neighborhood_cardinality_[i];
    }
  }
}

void NematicAlignmentSystem::EvaluateInteractionsWithVerletNeighborList(std::vector<Real> &system_state,
                                                                        const std::vector<Real> &k_prev,
                                                                        std::vector<Real> &k_next,
                                                                        Real k_coef,
                                                                        Real dt)
{
  std::fill(alignment_force_.begin(), alignment_force_.end(), 0.0);
  std::fill(neighborhood_cardinality_.begin(), neighborhood_cardinality_.end(), 0.0);

  Real x_i = 0.0, y_i = 0.0, phi_i = 0.0;
  Real x_j = 0.0, y_j = 0.0, phi_j = 0.0;
  Real dx = 0.0, dy = 0.0, dr_squared = 0.0;
  int i_cell = 0, i_cell_x = 0, i_cell_y = 0;
  int j_cell = 0, j_cell_x = 0, j_cell_y = 0;
  int ii = 0, j = 0, jj = 0;

  //if it is time to update the linked list and the Verlet list
  if (should_update_lists_)
  {
    //update the linked list
    for (int i = 0; i < linked_list_.size(); ++i)
    {
      linked_list_[i].clear();
    }

    //construct a linked list
#if defined(MPI_FAST_INTERACTION)
    for (int i : thread_.GetLoopIndices())
        {
            ii = kS * i;

            x_i   = system_state[ii    ] + k_coef * dt * k_prev[ii    ];
            y_i   = system_state[ii + 1] + k_coef * dt * k_prev[ii + 1];
//			phi_i = system_state[ii + 2] + k_coef * dt * k_prev[ii + 2];

            if (x_i < 0.0 || x_i >= x_size_ || y_i < 0.0 || y_i >= y_size_)
            {
                pbc_config_.ApplyPeriodicBoundaryConditions(x_i, y_i, x_i, y_i);
            }

            i_cell_x = int(x_i / x_size_ * num_subcells_x_);
            i_cell_y = int(y_i / y_size_ * num_subcells_y_);
            i_cell = i_cell_y * num_subcells_x_ + i_cell_x;

            pre_linked_list_[i] = i_cell;
        }

        if (mpi_thread_rank != mpi_root_rank)
        {
            MPI_Gather(&pre_linked_list_[mpi_thread_rank * mpi_num_particles_per_thread], mpi_num_particles_per_thread, MPI_INT, &pre_linked_list_[0], mpi_num_particles_per_thread, MPI_INT, mpi_root_rank, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Gather(MPI_IN_PLACE, mpi_num_particles_per_thread, MPI_INT, &pre_linked_list_[mpi_root_rank * mpi_num_particles_per_thread], mpi_num_particles_per_thread, MPI_INT, mpi_root_rank, MPI_COMM_WORLD);
        }
        MPI_Bcast(&pre_linked_list_[0], (int)pre_linked_list_.size(), MPI_INT, mpi_root_rank, MPI_COMM_WORLD);

        for (int i = 0; i < kN; ++i)
        {
            linked_list_[pre_linked_list_[i]].push_back(i);
        }
#else
    for (int i : thread_.GetLoopIndices())
    {
      ii = kS * i;

      x_i = system_state[ii] + k_coef * dt * k_prev[ii];
      y_i = system_state[ii + 1] + k_coef * dt * k_prev[ii + 1];
//			phi_i = system_state[ii + 2] + k_coef * dt * k_prev[ii + 2];

      if (x_i < 0.0 || x_i >= x_size_ || y_i < 0.0 || y_i >= y_size_)
      {
        pbc_config_.ApplyPeriodicBoundaryConditions(x_i, y_i, x_i, y_i);
      }

      i_cell_x = int(x_i / x_size_ * num_subcells_x_);
      i_cell_y = int(y_i / y_size_ * num_subcells_y_);
      i_cell = i_cell_y * num_subcells_x_ + i_cell_x;

      linked_list_[i_cell].push_back(i);
    }
#endif

    //update the Verlet list
    for (int i = 0; i < verlet_list_.size(); ++i)
    {
      verlet_list_[i].clear();
    }

    for (int i : thread_.GetLoopIndices())
    {
      ii = kS * i;

      x_i = system_state[ii] + k_coef * dt * k_prev[ii];
      y_i = system_state[ii + 1] + k_coef * dt * k_prev[ii + 1];
      phi_i = system_state[ii + 2] + k_coef * dt * k_prev[ii + 2];

      if (x_i < 0.0 || x_i >= x_size_ || y_i < 0.0 || y_i >= y_size_)
      {
        pbc_config_.ApplyPeriodicBoundaryConditions(x_i, y_i, x_i, y_i);
      }

      i_cell_x = int(x_i / x_size_ * num_subcells_x_);
      i_cell_y = int(y_i / y_size_ * num_subcells_y_);
      i_cell = i_cell_y * num_subcells_x_ + i_cell_x;

      for (int neighboring_cell = 0; neighboring_cell < neighboring_cells_.size(); ++neighboring_cell)
      {
        j_cell_x = i_cell_x + neighboring_cells_[neighboring_cell][0];
        j_cell_y = i_cell_y + neighboring_cells_[neighboring_cell][1];
        AdjustNeighboringCellToPeriodicBoundaries(j_cell_x, j_cell_y);
        j_cell = j_cell_y * num_subcells_x_ + j_cell_x;

        for (std::vector<int>::iterator neighbr_iter = linked_list_[j_cell].begin();
             neighbr_iter != linked_list_[j_cell].end(); ++neighbr_iter)
        {
          j = *neighbr_iter;
          jj = kS * j;

          if (j != i)
          {
            x_j = system_state[jj] + k_coef * dt * k_prev[jj];
            y_j = system_state[jj + 1] + k_coef * dt * k_prev[jj + 1];
            phi_j = system_state[jj + 2] + k_coef * dt * k_prev[jj + 2];

            if (x_j < 0.0 || x_j >= x_size_ || y_j < 0.0 || y_j >= y_size_)
            {
              pbc_config_.ApplyPeriodicBoundaryConditions(x_j, y_j, x_j, y_j);
            }

            pbc_config_.ClassAEffectiveParticleDistance(x_i, y_i, x_j, y_j, dx, dy);
            dr_squared = dx * dx + dy * dy;

            if (dr_squared <= r_max_ * r_max_)
            {
#if defined(MPI_FAST_INTERACTION)
              verlet_list_[i - mpi_thread_rank * mpi_num_particles_per_thread].push_back(j);
#else
              verlet_list_[i].push_back(j);
#endif
            }
          }
        } // neighbr_iter
      } // neighboring_cell
    } // i

    //reset the counter of the passed time steps
    should_update_lists_ = false;
  }

  //loop using the Verlet list
  for (int i : thread_.GetLoopIndices())
  {
    ii = kS * i;

    x_i = system_state[ii] + k_coef * dt * k_prev[ii];
    y_i = system_state[ii + 1] + k_coef * dt * k_prev[ii + 1];
    phi_i = system_state[ii + 2] + k_coef * dt * k_prev[ii + 2];

    if (x_i < 0.0 || x_i >= x_size_ || y_i < 0.0 || y_i >= y_size_)
    {
      pbc_config_.ApplyPeriodicBoundaryConditions(x_i, y_i, x_i, y_i);
    }

#if defined(MPI_FAST_INTERACTION)
    for (std::vector<int>::iterator neighbr_iter = verlet_list_[i - mpi_thread_rank * mpi_num_particles_per_thread].begin(); neighbr_iter != verlet_list_[i - mpi_thread_rank * mpi_num_particles_per_thread].end(); ++neighbr_iter)
#else
    for (std::vector<int>::iterator neighbr_iter = verlet_list_[i].begin(); neighbr_iter != verlet_list_[i].end();
         ++neighbr_iter)
#endif
    {
      j = *neighbr_iter;
      jj = kS * j;

      x_j = system_state[jj] + k_coef * dt * k_prev[jj];
      y_j = system_state[jj + 1] + k_coef * dt * k_prev[jj + 1];
      phi_j = system_state[jj + 2] + k_coef * dt * k_prev[jj + 2];

      if (x_j < 0.0 || x_j >= x_size_ || y_j < 0.0 || y_j >= y_size_)
      {
        pbc_config_.ApplyPeriodicBoundaryConditions(x_j, y_j, x_j, y_j);
      }

      pbc_config_.ClassAEffectiveParticleDistance(x_i, y_i, x_j, y_j, dx, dy);
      dr_squared = dx * dx + dy * dy;

      if (dr_squared <= r_min_ * r_min_)
      {
        alignment_force_[i] += std::sin(2.0 * (phi_j - phi_i));
        ++neighborhood_cardinality_[i];
      }
    } // neighbr_iter

    k_next[ii] = std::cos(phi_i);
    k_next[ii + 1] = std::sin(phi_i);
    if (neighborhood_cardinality_[i] != 0)
    {
      k_next[ii + 2] = sigma_ * alignment_force_[i] / neighborhood_cardinality_[i];
    }
  } // i

  if (last_coefficient_)
  {
//	Real max = 0.0;
//	Real dist = 0.0;
//	for (int i = 0; i < kN; ++i)
//	{
//		dist = 1.0 * dt;
//		if (max < dist)
//		{
//			max = dist;
//		}
//	}
    accumulated_displacement_ += 1.0 * dt;
//	if (mpi_thread_rank == mpi_root_rank)
//		std::cout << "accumulated_displacement_ = " << accumulated_displacement_ << ", r_max_ = " << r_max_ << ", r_min_ = " << r_min_ << std::endl;
    if (accumulated_displacement_ > 0.5 * (r_max_ - r_min_))
    {
      accumulated_displacement_ = 0.0;
      should_update_lists_ = true;
    }
  }
}

void NematicAlignmentSystem::AdjustNeighboringCellToPeriodicBoundaries(int &cell_x, int &cell_y)
{
  if (cell_x < 0)
  {
    cell_x = num_subcells_x_ - 1;
  } else if (cell_x >= num_subcells_x_)
  {
    cell_x = 0;
  }

  if (cell_y < 0)
  {
    cell_y = num_subcells_y_ - 1;
  } else if (cell_y >= num_subcells_y_)
  {
    cell_y = 0;
  }
}