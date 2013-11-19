//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// This class
#include "nitridation_simulation.h"

// GRINS
#include "grins/simulation_builder.h"
#include "grins/bc_handling_base.h"
#include "grins/reacting_low_mach_navier_stokes_bc_handling.h"
#include "grins/catalytic_wall_base.h"
#include "grins/antioch_chemistry.h"

// Antioch
#include "antioch/chemical_mixture.h"

namespace NitridationCalibration
{
  NitridationSimulation::NitridationSimulation( const GetPot& input,
					        GRINS::SimulationBuilder& sim_builder )
    : GRINS::Simulation(input,sim_builder)
  {
    return;
  }

  NitridationSimulation::~NitridationSimulation()
  {
    return;
  }

  void NitridationSimulation::set_gamma_CN( const double gamma )
  {
    std::tr1::shared_ptr<GRINS::Physics> physics = this->_multiphysics_system->get_physics( GRINS::reacting_low_mach_navier_stokes );

    GRINS::BCHandlingBase* bc_handler_base = physics->get_bc_handler();

    GRINS::ReactingLowMachNavierStokesBCHandling<GRINS::AntiochChemistry>* bc_handler =
      libmesh_cast_ptr<GRINS::ReactingLowMachNavierStokesBCHandling<GRINS::AntiochChemistry>*>(bc_handler_base);
    
    /*! \todo Need to generalize to more than 1 bc */
    
    const unsigned int bc_id = 3;

    GRINS::CatalyticWallBase<GRINS::AntiochChemistry>* func = bc_handler->get_catalytic_wall( bc_id );

    std::vector<libMesh::Real> params;
    params.push_back( gamma );

    func->set_catalycity_params( params );

    return;
  }

  void NitridationSimulation::reset_initial_guess( const NumericVector<libMesh::Real>& solution )
  {
    (*((this->_multiphysics_system)->solution)) = solution;
    (*((this->_multiphysics_system)->current_local_solution)) = solution;


    (*((this->_multiphysics_system)->solution)).close();
    (*((this->_multiphysics_system)->current_local_solution)).close();

    return;
  }

  void NitridationSimulation::solve()
  {
    this->_multiphysics_system->solve();

    return;
  }

  double NitridationSimulation::computed_mass_loss()
  {
    _multiphysics_system->assemble_qoi( libMesh::QoISet( *_multiphysics_system ) );
    return this->get_qoi_value(0);
  }

} // end namespace NitridationCalibration
