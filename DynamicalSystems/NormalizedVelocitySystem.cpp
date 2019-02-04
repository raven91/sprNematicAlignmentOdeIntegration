//
// Created by Nikita Kruk on 21.11.18.
//

#include "NormalizedVelocitySystem.hpp"

#include <eigen3/Eigen/Dense>

const Real rho = 0.1;
const Real U_0 = 0.1;

NormalizedVelocitySystem::NormalizedVelocitySystem(Thread &thread, PeriodicBoundaryConditionsConfiguration &pbc_config)
    :
    thread_(thread),
    pbc_config_(pbc_config),
    alignment_force_(4, std::vector<Real>(kN, 0.0)),
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
}

NormalizedVelocitySystem::~NormalizedVelocitySystem()
{
  alignment_force_.clear();
  neighborhood_cardinality_.clear();
  neighboring_cells_.clear();
  for (int i = 0; i < linked_list_.size(); ++i)
  {
    linked_list_[i].clear();
  }
  linked_list_.clear();
}

void NormalizedVelocitySystem::EvaluateRhs(std::vector<Real> &system_state,
                                           const std::vector<Real> &k_prev,
                                           std::vector<Real> &k_next,
                                           Real k_coef,
                                           Real dt)
{
//  EvaluateInteractionsWithLinkedList(system_state, k_prev, k_next, k_coef, dt);
  EvaluateInteractionsWithAllPairs(system_state, k_prev, k_next, k_coef, dt);
}

void NormalizedVelocitySystem::EvaluateInteractionsWithAllPairs(std::vector<Real> &system_state,
                                                                const std::vector<Real> &k_prev,
                                                                std::vector<Real> &k_next,
                                                                Real k_coef,
                                                                Real dt)
{
  std::fill(alignment_force_[0].begin(), alignment_force_[0].end(), 0.0);
  std::fill(alignment_force_[1].begin(), alignment_force_[1].end(), 0.0);
  std::fill(alignment_force_[2].begin(), alignment_force_[2].end(), 0.0);
  std::fill(alignment_force_[3].begin(), alignment_force_[3].end(), 0.0);
  std::fill(neighborhood_cardinality_.begin(), neighborhood_cardinality_.end(), 0.0);

  Real x_i = 0.0, y_i = 0.0, u_x_i = 0.0, u_y_i = 0.0;
  Real x_j = 0.0, y_j = 0.0, u_x_j = 0.0, u_y_j = 0.0;
  Real dx = 0.0, dy = 0.0, dr_squared = 0.0;

  for (int i : thread_.GetLoopIndices())
  {
    int ii = kS * i;

    x_i = system_state[ii] + k_coef * dt * k_prev[ii];
    y_i = system_state[ii + 1] + k_coef * dt * k_prev[ii + 1];
    u_x_i = system_state[ii + 2] + k_coef * dt * k_prev[ii + 2];
    u_y_i = system_state[ii + 3] + k_coef * dt * k_prev[ii + 3];

    //j ~ i
    for (int j = 0; j < kN; ++j)
    {
      int jj = kS * j;

      if (j != i)
      {
        x_j = system_state[jj] + k_coef * dt * k_prev[jj];
        y_j = system_state[jj + 1] + k_coef * dt * k_prev[jj + 1];
        u_x_j = system_state[jj + 2] + k_coef * dt * k_prev[jj + 2];
        u_y_j = system_state[jj + 3] + k_coef * dt * k_prev[jj + 3];

        pbc_config_.ClassAEffectiveParticleDistance(x_i, y_i, x_j, y_j, dx, dy);
//        dr_squared = dx * dx + dy * dy;

//        Real rho_squared = rho * rho;
//        if (dr_squared <= rho_squared)
        {
//          Eigen::Vector4d dr_0 = YukawaPotentialForce(x_i, y_i, x_j, y_j, u_x_j, u_y_j);
          // TODO: sign is determined by which variable is differentiated
//          Eigen::Vector4d dr_1 = -YukawaPotentialForce(x_j, y_j, x_i, y_i, u_x_i, u_y_i);
          Eigen::Vector2d dr_0 = YukawaForce(x_i, y_i, x_j, y_j, u_x_i, u_y_i);
          Eigen::Vector2d dr_1 = YukawaForce(x_i, y_i, x_j, y_j, u_x_j, u_y_j);

          alignment_force_[0][i] += 0.5 * (dr_0(0) + dr_1(0));
          alignment_force_[1][i] += 0.5 * (dr_0(1) + dr_1(1));
//          alignment_force_[2][i] += 0.5 * (dr_0(2) + dr_1(2));
//          alignment_force_[3][i] += 0.5 * (dr_0(3) + dr_1(3));
          ++neighborhood_cardinality_[i];
        }
      }
    } // j

    k_next[ii] = 1.0 * u_x_i;
    k_next[ii + 1] = 1.0 * u_y_i;
    if (neighborhood_cardinality_[i] != 0.0)
    {
//      k_next[ii] = 1.0 * u_x_i + alignment_force_[0][i] / neighborhood_cardinality_[i];
//      k_next[ii + 1] = 1.0 * u_y_i + alignment_force_[1][i] / neighborhood_cardinality_[i];
      Real alpha = 1.0, beta = 1.0;
      Real norm_squared = u_x_i * u_x_i + u_y_i * u_y_i;
      k_next[ii + 2] =
          (alpha - beta * norm_squared) * u_x_i + alignment_force_[0][i] / neighborhood_cardinality_[i];
      k_next[ii + 3] =
          (alpha - beta * norm_squared) * u_y_i + alignment_force_[1][i] / neighborhood_cardinality_[i];
    }
  } // i
}

void NormalizedVelocitySystem::EvaluateInteractionsWithLinkedList(std::vector<Real> &system_state,
                                                                  const std::vector<Real> &k_prev,
                                                                  std::vector<Real> &k_next,
                                                                  Real k_coef,
                                                                  Real dt)
{
  std::fill(alignment_force_[0].begin(), alignment_force_[0].end(), 0.0);
  std::fill(alignment_force_[1].begin(), alignment_force_[1].end(), 0.0);
  std::fill(neighborhood_cardinality_.begin(), neighborhood_cardinality_.end(), 0.0);

  Real x_i = 0.0, y_i = 0.0, u_x_i = 0.0, u_y_i = 0.0;
  Real x_j = 0.0, y_j = 0.0, u_x_j = 0.0, u_y_j = 0.0;
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
    u_x_i = system_state[ii + 2] + k_coef * dt * k_prev[ii + 2];
    u_y_i = system_state[ii + 3] + k_coef * dt * k_prev[ii + 3];

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
          u_x_j = system_state[jj + 2] + k_coef * dt * k_prev[jj + 2];
          u_y_j = system_state[jj + 3] + k_coef * dt * k_prev[jj + 3];

          if (x_j < 0.0 || x_j >= x_size_ || y_j < 0.0 || y_j >= y_size_)
          {
            pbc_config_.ApplyPeriodicBoundaryConditions(x_j, y_j, x_j, y_j);
          }

          pbc_config_.ClassAEffectiveParticleDistance(x_i, y_i, x_j, y_j, dx, dy);
          dr_squared = dx * dx + dy * dy;

          Real rho_squared = rho * rho;
//        if (dr_squared <= rho_squared)
          {
            Eigen::Vector4d dr_0 = YukawaPotentialForce(x_i, y_i, x_j, y_j, u_x_j, u_y_j);
            // TODO: sign is determined by which variable is differentiated
            Eigen::Vector4d dr_1 = -YukawaPotentialForce(x_j, y_j, x_i, y_i, u_x_i, u_y_i);

            alignment_force_[0][i] += 0.5 * (dr_0(0) + dr_1(0));
            alignment_force_[1][i] += 0.5 * (dr_0(1) + dr_1(1));
            alignment_force_[2][i] += 0.5 * (dr_0(2) + dr_1(2));
            alignment_force_[3][i] += 0.5 * (dr_0(3) + dr_1(3));
            ++neighborhood_cardinality_[i];
          }
        }
      } // neighbr_iter
    } // neighboring_cell

    k_next[ii] = 1.0 * u_x_i;
    k_next[ii + 1] = 1.0 * u_y_i;
    if (neighborhood_cardinality_[i] != 0.0)
    {
      k_next[ii + 2] = 1.0 * alignment_force_[0][i] / neighborhood_cardinality_[i];
      k_next[ii + 3] = 1.0 * alignment_force_[1][i] / neighborhood_cardinality_[i];
    }
  } // i
}

void NormalizedVelocitySystem::AdjustNeighboringCellToPeriodicBoundaries(int &cell_x, int &cell_y)
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

Eigen::Vector4d NormalizedVelocitySystem::YukawaPotentialForce(Real x_i,
                                                               Real y_i,
                                                               Real x_j,
                                                               Real y_j,
                                                               Real u_x,
                                                               Real u_y)
{
  Real lambda_x = 1.0, lambda_y = 1.0, lambda = 1.0;
  Eigen::Vector2d dr(x_j - x_i, y_j - y_i);
  Eigen::Matrix2d rotation;
  rotation << u_y, u_x, -u_x, u_y;
  Eigen::Matrix2d inverse_rotation = rotation.inverse().eval();
  Eigen::Vector2d dr_rotated = inverse_rotation * dr;

  Real s = std::sqrt(
      dr_rotated(0) * dr_rotated(0) / (lambda_x * lambda_x) + dr_rotated(1) * dr_rotated(1) / (lambda_y * lambda_y));
  Real s_recip = 1.0 / s;
  Real strength = -s_recip * (s_recip + 1.0 / lambda) * std::exp(-s / lambda);
  Eigen::Vector4d gradient;
  gradient(0) = -s_recip * (dr_rotated(0) * u_y / (lambda_x * lambda_x) - dr_rotated(1) * u_x / (lambda_y * lambda_y));
  gradient(1) = -s_recip * (dr_rotated(0) * u_x / (lambda_x * lambda_x) + dr_rotated(1) * u_y / (lambda_y * lambda_y));
  gradient(2) =
      -s_recip * (dr_rotated(0) * dr(1) / (lambda_x * lambda_x) - dr_rotated(1) * dr(0) / (lambda_y * lambda_y));
  gradient(3) =
      -s_recip * (dr_rotated(0) * dr(0) / (lambda_x * lambda_x) + dr_rotated(1) * dr(1) / (lambda_y * lambda_y));

  gradient.head(2) = -U_0 * rotation * (strength * gradient.head(2));
  gradient.tail(2) = -U_0 * rotation * (strength * gradient.tail(2));
  return gradient;
}

Eigen::Vector2d NormalizedVelocitySystem::YukawaForce(Real x_i, Real y_i, Real x_j, Real y_j, Real u_x, Real u_y)
{
  Real lambda_x = 1.0, lambda_y = 4.0, lambda = 1.0;
  Eigen::Vector2d dr(x_j - x_i, y_j - y_i);
  Eigen::Matrix2d rotation;
  rotation << u_y, u_x, -u_x, u_y;
  Eigen::Matrix2d inverse_rotation = rotation.inverse().eval();
  Eigen::Vector2d dr_rotated = inverse_rotation * dr;

  Real s = std::sqrt(
      dr_rotated(0) * dr_rotated(0) / (lambda_x * lambda_x) + dr_rotated(1) * dr_rotated(1) / (lambda_y * lambda_y));
  Real s_recip = 1.0 / s;
  Real force = U_0 * s_recip * (s_recip + 1.0 / lambda) * std::exp(-s / lambda);

  return (-force * dr.normalized());
}