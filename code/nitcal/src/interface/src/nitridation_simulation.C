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
#include "grins/catalytic_wall.h"

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

    GRINS::BCHandlingBase* bc_handler = physics->get_bc_handler();

    // CN part
    {
      const unsigned int bc_id = 3;
      const GRINS::VariableIndex var_id = _multiphysics_system->variable_number("w_CN");

      std::tr1::shared_ptr< GRINS::NeumannFuncObj > raw_func =
	bc_handler->get_neumann_bound_func( bc_id, var_id );

      GRINS::CatalyticWall<Antioch::ChemicalMixture<libMesh::Real> >* func = libmesh_cast_ptr<GRINS::CatalyticWall<Antioch::ChemicalMixture<libMesh::Real> >* >( raw_func.get() );

      func->set_gamma( gamma );
    }

    // N part
    {
      const unsigned int bc_id = 3;
      const GRINS::VariableIndex var_id = this->_multiphysics_system->variable_number("w_N");

      std::tr1::shared_ptr< GRINS::NeumannFuncObj > raw_func =
	bc_handler->get_neumann_bound_func( bc_id, var_id );

      GRINS::CatalyticWall<Antioch::ChemicalMixture<libMesh::Real> >* func = libmesh_cast_ptr<GRINS::CatalyticWall<Antioch::ChemicalMixture<libMesh::Real> >* >( raw_func.get() );

      func->set_gamma( -gamma );
    }

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
    return this->get_qoi(0);
  }

} // end namespace NitridationCalibration
