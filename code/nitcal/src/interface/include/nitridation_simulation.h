//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef NITCAL_NITRIDATION_SIMULATION_H
#define NITCAL_NITRIDATION_SIMULATION_H

// GRINS
#include "grins/simulation.h"

namespace NitridationCalibration
{

  class NitridationSimulation : public GRINS::Simulation
  {
  public:

    NitridationSimulation( const GetPot& input,
			   GRINS::SimulationBuilder& sim_builder );

    virtual ~NitridationSimulation();

    void set_gamma_CN_params( const std::vector<double>& gamma_CN_params );

    void set_gamma_N_params( const std::vector<double>& gamma_N_params );

    void reset_initial_guess( const NumericVector<libMesh::Real>& solution );

    void solve();

    double computed_mass_loss();

    double computed_average_n();

  protected:

    unsigned int _mass_loss_qoi_index;

    unsigned int _mass_loss_catalytic_qoi_index;

    unsigned int _average_n_qoi_index;

  };

} // end namespace NitridationCalibration

#endif // NITCAL_NITRIDATION_SIMULATION_H
