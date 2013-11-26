//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef NITCAL_SIMULATION_INTERFACE_H
#define NITCAL_SIMULATION_INTERFACE_H

// Nitridation
#include "nitridation_simulation.h"

namespace NitridationCalibration
{

  class SimulationInterface
  {
  public:

    SimulationInterface( int argc, char** argv, MPI_Comm mpi_comm,
			 const GetPot& input );
    ~SimulationInterface();

    void set_gamma_CN_params( const std::vector<double>& gamma_CN_params ) const;

    void set_gamma_N_params( const std::vector<double>& gamma_N_params ) const;

    void reset_initial_guess() const;

    void solve() const;
    
    double computed_mass_loss() const;

    double computed_average_n() const;

  private:

    libMesh::LibMeshInit _libmesh_init;

    NitridationSimulation* _simulation;

    libMesh::AutoPtr<libMesh::NumericVector<libMesh::Real> > _cached_initial_guess;

  };

} // namespace Damage

#endif // ONE_D_DAMAGE_INTERFACE_H
