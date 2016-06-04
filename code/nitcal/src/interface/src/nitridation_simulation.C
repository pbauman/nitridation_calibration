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
#include "grins/catalytic_wall_base.h"
#include "grins/antioch_chemistry.h"
#include "grins/neumann_bc_container.h"
#include "grins/gas_solid_catalytic_wall.h"
#include "grins/gas_recombination_catalytic_wall.h"
#include "grins/antioch_chemistry.h"

// Antioch
#include "antioch/chemical_mixture.h"

namespace NitridationCalibration
{
  NitridationSimulation::NitridationSimulation( const GetPot& input,
                                                GetPot& command_line,
					        GRINS::SimulationBuilder& sim_builder,
                                                const libMesh::Parallel::Communicator& comm )
    : GRINS::Simulation(input,command_line,sim_builder,comm),
      _gas_recomb_idx(std::numeric_limits<unsigned int>::max()),
      _gas_solid_idx(std::numeric_limits<unsigned int>::max())
  {
    /* Search for the QoI's and cache their indices for later use */
    // First get the DifferentiableQoI and cast to CompositeQoI
    libMesh::DifferentiableQoI* qoi_base = _multiphysics_system->get_qoi();
    GRINS::CompositeQoI* qois = libMesh::libmesh_cast_ptr<GRINS::CompositeQoI*>( qoi_base );

    const unsigned int n_qois = qois->n_qois();

    for( unsigned int q = 0; q < n_qois; q++ )
      {
        const GRINS::QoIBase& qoi = qois->get_qoi(q);
        const std::string& qoi_name = qoi.name();

        if( qoi_name == average_N_mole_fraction )
          _average_n_qoi_index = q;

        else if( qoi_name == mass_loss_catalytic )
          _mass_loss_catalytic_qoi_index = q;

        else
          libmesh_error_msg("Error: Invalid qoi_name " << qoi_name << std::endl);
      }

    // Get the Neumann BCs and search for the boundary id corresponding to
    // each of the catalytic walls and cache the index later for resetting
    // the parameters
    const std::vector<GRINS::SharedPtr<GRINS::NeumannBCContainer> > & neumann_bcs =
      this->_multiphysics_system->get_neumann_bcs();

    for( unsigned int i = 0; i < neumann_bcs.size(); i++ )
      {
        if( neumann_bcs[i]->has_bc_id(2) )
          _gas_recomb_idx = i;

        if( neumann_bcs[i]->has_bc_id(3) )
          _gas_solid_idx = i;
      }

    if( _gas_recomb_idx == std::numeric_limits<unsigned int>::max() )
      libmesh_error_msg("ERROR: Could not find idx for GasRecomb BC!");

    if( _gas_solid_idx == std::numeric_limits<unsigned int>::max() )
      libmesh_error_msg("ERROR: Could not find idx for GasSolid BC!");
  }

  void NitridationSimulation::set_gamma_CN_params( const std::vector<double>& gamma_CN_params )
  {
    std::vector<GRINS::SharedPtr<GRINS::NeumannBCContainer> > & neumann_bcs =
      this->_multiphysics_system->get_neumann_bcs();

    GRINS::SharedPtr<GRINS::NeumannBCAbstract> bc_base =
      neumann_bcs[_gas_solid_idx]->get_func();

    GRINS::GasSolidCatalyticWall<GRINS::AntiochChemistry>* wall =
      libMesh::cast_ptr<GRINS::GasSolidCatalyticWall<GRINS::AntiochChemistry>*>( bc_base.get() );

    wall->set_catalycity_params( gamma_CN_params );
  }

  void NitridationSimulation::set_gamma_N_params( const std::vector<double>& gamma_N_params )
  {
    std::vector<GRINS::SharedPtr<GRINS::NeumannBCContainer> > & neumann_bcs =
      this->_multiphysics_system->get_neumann_bcs();

    GRINS::SharedPtr<GRINS::NeumannBCAbstract> bc_base =
      neumann_bcs[_gas_recomb_idx]->get_func();

    GRINS::GasRecombinationCatalyticWall<GRINS::AntiochChemistry>* wall =
      libMesh::cast_ptr<GRINS::GasRecombinationCatalyticWall<GRINS::AntiochChemistry>*>( bc_base.get() );

    wall->set_catalycity_params( gamma_N_params );
  }

  void NitridationSimulation::reset_initial_guess( const libMesh::NumericVector<libMesh::Real>& solution )
  {
    (*((this->_multiphysics_system)->solution)) = solution;
    (*((this->_multiphysics_system)->current_local_solution)) = solution;


    (*((this->_multiphysics_system)->solution)).close();
    (*((this->_multiphysics_system)->current_local_solution)).close();
  }

  void NitridationSimulation::solve()
  {
    this->_multiphysics_system->solve();
  }

  double NitridationSimulation::computed_mass_loss()
  {
    _multiphysics_system->assemble_qoi( libMesh::QoISet( *_multiphysics_system ) );

    /* Absolute value because this will be a negative quantity, but the data
       to which we are comparing are all given as postive values */
    return std::fabs(this->get_qoi_value(_mass_loss_catalytic_qoi_index));
  }

  double NitridationSimulation::computed_average_n()
  {
    _multiphysics_system->assemble_qoi( libMesh::QoISet( *_multiphysics_system ) );
    return this->get_qoi_value(_average_n_qoi_index);
  }

} // end namespace NitridationCalibration
