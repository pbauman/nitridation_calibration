//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// This class
#include "nitridation_simulation.h"

// NitCal
#include "qoi_names.h"

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
    /* Search for the QoI's and cache their indices for later use */
    // First get the DifferentiableQoI and cast to CompositeQoI
    libMesh::DifferentiableQoI* qoi_base = _multiphysics_system->get_qoi();
    GRINS::CompositeQoI* qois = libmesh_cast_ptr<GRINS::CompositeQoI*>( qoi_base );

    const unsigned int n_qois = qois->n_qois();

    for( unsigned int q = 0; q < n_qois; q++ )
      {
        const GRINS::QoIBase& qoi = qois->get_qoi(q);

        const std::string& qoi_name = qoi.name();

        if( qoi_name == average_N_mole_fraction )
          {
            _average_n_qoi_index = q;
          }
        else if( qoi_name == mass_loss_catalytic )
          {
            _mass_loss_catalytic_qoi_index = q;
          }
        else
          {
            std::cerr << "Error: Invalid qoi_name " << qoi_name << std::endl;
            libmesh_error();
          }
      }
    
    return;
  }

  NitridationSimulation::~NitridationSimulation()
  {
    return;
  }

  void NitridationSimulation::set_gamma_CN_params( const std::vector<double>& gamma_CN_params )
  {
    std::tr1::shared_ptr<GRINS::Physics> physics = this->_multiphysics_system->get_physics( GRINS::reacting_low_mach_navier_stokes );

    GRINS::BCHandlingBase* bc_handler_base = physics->get_bc_handler();

    GRINS::ReactingLowMachNavierStokesBCHandling<GRINS::AntiochChemistry>* bc_handler =
      libmesh_cast_ptr<GRINS::ReactingLowMachNavierStokesBCHandling<GRINS::AntiochChemistry>*>(bc_handler_base);
    
    /*! \todo Need to generalize to more than 1 bc */
    
    const unsigned int bc_id = 3;

    GRINS::CatalyticWallBase<GRINS::AntiochChemistry>* func = bc_handler->get_catalytic_wall( bc_id );

    func->set_catalycity_params( gamma_CN_params );

    return;
  }

  void NitridationSimulation::set_gamma_N_params( const std::vector<double>& gamma_N_params )
  {
    std::tr1::shared_ptr<GRINS::Physics> physics = this->_multiphysics_system->get_physics( GRINS::reacting_low_mach_navier_stokes );

    GRINS::BCHandlingBase* bc_handler_base = physics->get_bc_handler();

    GRINS::ReactingLowMachNavierStokesBCHandling<GRINS::AntiochChemistry>* bc_handler =
      libmesh_cast_ptr<GRINS::ReactingLowMachNavierStokesBCHandling<GRINS::AntiochChemistry>*>(bc_handler_base);
    
    /*! \todo Need to generalize to more than 1 bc */
    
    const unsigned int bc_id = 2;

    GRINS::CatalyticWallBase<GRINS::AntiochChemistry>* func = bc_handler->get_catalytic_wall( bc_id );

    func->set_catalycity_params( gamma_N_params );

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
    return this->get_qoi_value(_mass_loss_catalytic_qoi_index);
  }

  double NitridationSimulation::computed_average_n()
  {
    _multiphysics_system->assemble_qoi( libMesh::QoISet( *_multiphysics_system ) );
    return this->get_qoi_value(_average_n_qoi_index);
  }

} // end namespace NitridationCalibration
