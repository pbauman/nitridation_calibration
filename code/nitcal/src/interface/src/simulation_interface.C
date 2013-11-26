//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// This class
#include "simulation_interface.h"

// NitCal
#include "nitcal_bc_factory.h"
#include "qoi_factory.h"

// GRINS
#include "grins/simulation_builder.h"
#include "grins/bc_factory.h"
#include "grins/qoi_factory.h"

namespace NitridationCalibration
{
  // Function for getting initial temperature field
  libMesh::Real initial_values( const libMesh::Point& p,
				const libMesh::Parameters &params, 
				const std::string& system_name,
				const std::string& unknown_name );

  SimulationInterface::SimulationInterface( int argc, char** argv,
					    MPI_Comm mpi_comm,
					    const GetPot& input )
    : _libmesh_init(argc,argv,mpi_comm),
      _cached_initial_guess( libMesh::NumericVector<libMesh::Real>::build() )
  {
    GRINS::SimulationBuilder sim_builder;
    
    std::tr1::shared_ptr<GRINS::BoundaryConditionsFactory> bc_factory( new NitridationCalibration::BoundaryConditionsFactory(input) );

    sim_builder.attach_bc_factory( bc_factory );
    
    std::tr1::shared_ptr<GRINS::QoIFactory> qoi_factory( new NitridationCalibration::QoIFactory );
    
    sim_builder.attach_qoi_factory( qoi_factory );

    _simulation = new NitridationSimulation( input, sim_builder );

    // Project initial solution
    std::string restart_file = input( "restart-options/restart_file", "none" );

    std::tr1::shared_ptr<NitridationCalibration::TubeTempBC> wall_temp;
  
    if( restart_file == "none" )
      {
	// Asssign initial temperature value
	std::string system_name = input( "screen-options/system_name", "GRINS" );
	std::tr1::shared_ptr<libMesh::EquationSystems> es = _simulation->get_equation_system();
	const libMesh::System& system = es->get_system(system_name);
      
	Parameters &params = es->parameters;

	Real& w_N2 = params.set<Real>( "w_N2" );
	w_N2 = input( "Physics/ReactingLowMachNavierStokes/bound_species_1", 0.0, 0 );
      
	Real& w_N = params.set<Real>( "w_N" );
	w_N = input( "Physics/ReactingLowMachNavierStokes/bound_species_1", 0.0, 1 );

	wall_temp.reset( new NitridationCalibration::TubeTempBC( input ) );
	std::tr1::shared_ptr<NitridationCalibration::TubeTempBC>& dummy = params.set<std::tr1::shared_ptr<NitridationCalibration::TubeTempBC> >( "wall_temp" );
	dummy = wall_temp;

	system.project_solution( initial_values, NULL, params );

	_cached_initial_guess->init( *(system.solution.get()), true );

	*(_cached_initial_guess.get()) = *(system.solution.get());
      }
    
    return;
  }

  SimulationInterface::~SimulationInterface()
  {
    delete _simulation;
    return;
  }

  void SimulationInterface::set_gamma_CN_params( const std::vector<double>& gamma_CN_params ) const
  {
    _simulation->set_gamma_CN_params(gamma_CN_params);
    return;
  }

  void SimulationInterface::set_gamma_N_params( const std::vector<double>& gamma_N_params ) const
  {
    _simulation->set_gamma_N_params(gamma_N_params);
    return;
  }

  void SimulationInterface::reset_initial_guess() const
  {
    _simulation->reset_initial_guess( *_cached_initial_guess );
    return;
  }

  void SimulationInterface::solve() const
  {
    _simulation->solve();
    return;
  }
  
  double SimulationInterface::computed_mass_loss() const
  {
    return _simulation->computed_mass_loss();
  }

  double SimulationInterface::computed_average_n() const
  {
    return _simulation->computed_average_n();
  }

  libMesh::Real initial_values( const libMesh::Point& p,
				const libMesh::Parameters &params, 
				const std::string& ,
				const std::string& unknown_name )
  {
    libMesh::Real value = 0.0;

    if( unknown_name == "w_N2" )
      value = params.get<libMesh::Real>("w_N2");

    else if( unknown_name == "w_N" )
      value = params.get<libMesh::Real>("w_N");

    else if( unknown_name == "T" )
      {
	value = (*params.get<std::tr1::shared_ptr<NitridationCalibration::TubeTempBC> >( "wall_temp" ))(p);
	//value = 1200.0;
      }

    else if( unknown_name == "u" )
      //value = 40.0;
      value = 0.0;

    else
      value = 0.0;

    return value;
  }

} // end namespace NitridationCalibration
