//
// Created by Nikita Kruk on 16.11.18.
//

#include "BinaryObserver.hpp"

#include <sstream>
#include <cassert>
#include <cmath>
#include <complex>
#include <numeric> // std::accumulate
#include <iostream>
#include <cmath>

BinaryObserver::BinaryObserver(Thread &thread,
                               Real sigma,
                               Real rho,
                               Real D_T,
                               Real D_R,
                               PeriodicBoundaryConditionsConfiguration &pbc_config,
                               Real dt,
                               int trial) :
    thread_(thread),
    pbc_config_(pbc_config),
    output_time_counter_{0, 0},
    output_time_threshold_{20, 1} // mod 1 - save at every dt
{
  if (thread.IsRoot())
  {
    integration_step_timer_ = std::chrono::system_clock::now();

#if defined(__linux__) && defined(LICHTENBERG)
    std::string folder("/work/scratch/nk59zoce/cpp/sprNematicAlignmentOdeIntegration/");
#elif defined(__linux__) && defined(BCS_CLUSTER)
    std::string folder("/home/nkruk/cpp/sprNematicAlignmentOdeIntegration/output/");
#elif defined(__linux__)
    std::string folder("/home/nikita/Documents/sprNematicAlignmentOdeIntegration/");
#elif defined(__APPLE__)
    std::string folder("/Users/nikita/Documents/Projects/spr/sprNematicAlignmentOdeIntegration/");
#endif
    std::ostringstream simulation_file_name_buffer;
    simulation_file_name_buffer << folder << "sigma_" << sigma << "_rho_" << rho << "_D_T_" << D_T << "_D_R_" << D_R
                                << "_N_" << kN << "_" << thread_.GetRank() << "_" << trial << ".bin";
    simulation_file_name_ = simulation_file_name_buffer.str();
    std::remove(simulation_file_name_.c_str());
    simulation_file_.open(simulation_file_name_, std::ios::binary | std::ios::out | std::ios::app);
    assert(simulation_file_.is_open());

    std::ostringstream summary_statistics_file_name_buffer;
    summary_statistics_file_name_buffer << folder << "sigma_" << sigma << "_rho_" << rho << "_D_T_" << D_T << "_D_R_"
                                        << D_R << "_N_" << kN << "_" << thread_.GetRank() << "_" << trial << ".txt";
    summary_statistics_file_name_ = summary_statistics_file_name_buffer.str();
    std::remove(summary_statistics_file_name_.c_str());
    summary_statistics_file_.open(summary_statistics_file_name_, std::ios::out | std::ios::app);
    assert(summary_statistics_file_.is_open());
  }
}

BinaryObserver::~BinaryObserver()
{
  if (thread_.IsRoot())
  {
    if (simulation_file_.is_open())
    {
      simulation_file_.close();
    }
    if (summary_statistics_file_.is_open())
    {
      summary_statistics_file_.close();
    }
  }
}

void BinaryObserver::SaveSystemState(std::vector<Real> &system_state, Real t)
{
  if (thread_.IsRoot())
  {
#if defined(CARTESIAN_REPRESENTATION)
//    for (int i = 0; i < kN; ++i)
//    {
//      Real norm = std::sqrt(
//          system_state[kS * i + 2] * system_state[kS * i + 2] + system_state[kS * i + 3] * system_state[kS * i + 3]);
//      system_state[kS * i + 2] /= norm;
//      system_state[kS * i + 3] /= norm;
//    }
    std::vector<float> polar_system_state(3 * kN, 0.0);
    for (int i = 0; i < kN; ++i)
    {
      polar_system_state[3 * i] = system_state[kS * i];
      polar_system_state[3 * i + 1] = system_state[kS * i + 1];
      polar_system_state[3 * i + 2] = std::atan2(system_state[kS * i + 3], system_state[kS * i + 2]);
    } // i
#endif
    if (!(output_time_counter_[0] % output_time_threshold_[0]))
    {
      std::chrono::duration<Real> elapsed_seconds = std::chrono::system_clock::now() - integration_step_timer_;
      std::cout << "value at t = " << t << " integrated in " << elapsed_seconds.count() << "s" << std::endl;
      integration_step_timer_ = std::chrono::system_clock::now();

      float t_float = t;
      simulation_file_.write((char *) (&t_float), sizeof(float));
#if defined(CARTESIAN_REPRESENTATION)
      simulation_file_.write((char *) &polar_system_state[0], 3 * kN * sizeof(float));
#else
      std::vector<float> system_state_float(system_state.begin(), system_state.end());
      simulation_file_.write((char *) (&system_state_float[0]), kS * kN * sizeof(float));
#endif
    }
    ++output_time_counter_[0];
  }
}

void BinaryObserver::SaveSummaryStatistics(std::vector<Real> &system_state, Real t)
{
  if (thread_.IsRoot())
  {
    if (!(output_time_counter_[1] % output_time_threshold_[1]))
    {
      float t_float = t;
      std::complex<float> polar_order_parameter(0.0f, 0.0f);
      std::complex<float> nematic_order_parameter(0.0f, 0.0f);
      for (int i = 0; i < kN; ++i)
      {
#if defined(CARTESIAN_REPRESENTATION)
        polar_order_parameter +=
            std::complex<float>((float) system_state[kS * i + 2], (float) system_state[kS * i + 3]);
        nematic_order_parameter +=
            std::complex<float>(float(2.0 * system_state[kS * i + 2]), float(2.0 * system_state[kS * i + 3]));
#else
        polar_order_parameter +=
            std::complex<float>((float) std::cos(system_state[kS * i + 2]), (float) std::sin(system_state[kS * i + 2]));
        nematic_order_parameter += std::complex<float>((float) std::cos(2.0 * system_state[kS * i + 2]),
                                                       (float) std::sin(2.0 * system_state[kS * i + 2]));
#endif

      } // i
      polar_order_parameter /= kN;
      nematic_order_parameter /= kN;

      summary_statistics_file_ << t_float << '\t' << std::abs(polar_order_parameter) << '\t'
                               << std::arg(polar_order_parameter) << '\t'
                               << std::abs(nematic_order_parameter) << '\t'
                               << std::arg(nematic_order_parameter);
      summary_statistics_file_ << std::endl;
    }
    ++output_time_counter_[1];
  }
}