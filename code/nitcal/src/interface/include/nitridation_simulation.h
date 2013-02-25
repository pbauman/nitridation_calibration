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

    ~NitridationSimulation();

    void set_gamma_CN( const double gamma );

    void reset_initial_guess( const NumericVector<libMesh::Real>& solution );

    void solve();

    double computed_mass_loss();

  protected:

  };

} // end namespace NitridationCalibration

#endif // NITCAL_NITRIDATION_SIMULATION_H
