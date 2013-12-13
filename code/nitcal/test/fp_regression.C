//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

#include "grins_config.h"

#include <iostream>

// GRINS
#include "grins/simulation.h"
#include "grins/simulation_builder.h"

// libMesh
#include "libmesh/parallel.h"
#include "libmesh/exact_solution.h"

// NitCal
#include "nitcal_bc_factory.h"
#include "qoi_factory.h"

// Function for getting initial temperature field
Real initial_values( const Point& p, const Parameters &params, 
		     const std::string& system_name, const std::string& unknown_name );

int run( int argc, char* argv[], const GetPot& input );

int main(int argc, char* argv[])
{
  // Check command line count.
  if( argc < 3 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify libMesh input file and exact solution file." << std::endl;
      exit(1); // TODO: something more sophisticated for parallel runs?
    }

  // libMesh input file should be first argument
  std::string libMesh_input_filename = argv[1];
  
  // Create our GetPot object.
  GetPot libMesh_inputfile( libMesh_input_filename );

  int return_flag = run(argc,argv,libMesh_inputfile);

  return return_flag;
}

int run( int argc, char* argv[], const GetPot& input )
{
  // Initialize libMesh library.
  LibMeshInit libmesh_init(argc, argv);
 
  GRINS::SimulationBuilder sim_builder;

  std::tr1::shared_ptr<GRINS::BoundaryConditionsFactory> bc_factory( new NitridationCalibration::BoundaryConditionsFactory(input) );

  sim_builder.attach_bc_factory( bc_factory );

  std::tr1::shared_ptr<GRINS::QoIFactory> qoi_factory( new NitridationCalibration::QoIFactory );

  sim_builder.attach_qoi_factory( qoi_factory );

  GRINS::Simulation grins( input,
			   sim_builder );

  //FIXME: We need to move this to within the Simulation object somehow...
  std::string restart_file = input( "restart-options/restart_file", "none" );

  std::tr1::shared_ptr<NitridationCalibration::TubeTempBC> wall_temp;
  
  std::string system_name = input( "screen-options/system_name", "GRINS" );

  if( restart_file == "none" )
    {
      // Asssign initial temperature value
      std::tr1::shared_ptr<libMesh::EquationSystems> es = grins.get_equation_system();
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
    }

  grins.run();

  // Get equation systems to create ExactSolution object
  std::tr1::shared_ptr<EquationSystems> es = grins.get_equation_system();

  //es->write("foobar.xdr");

  // Create Exact solution object and attach exact solution quantities
  ExactSolution exact_sol(*es);
  
  EquationSystems es_ref( es->get_mesh() );

  // Filename of file where comparison solution is stashed
  std::string solution_file = std::string(argv[2]);
  es_ref.read( solution_file );

  exact_sol.attach_reference_solution( &es_ref );
  
  // Compute error and get it in various norms
  exact_sol.compute_error(system_name, "u");
  exact_sol.compute_error(system_name, "v");

  if( (es->get_mesh()).mesh_dimension() == 3 )
    exact_sol.compute_error(system_name, "w");

  exact_sol.compute_error(system_name, "p");
  exact_sol.compute_error(system_name, "T");
  exact_sol.compute_error(system_name, "w_N2");
  exact_sol.compute_error(system_name, "w_N");
  exact_sol.compute_error(system_name, "w_CN");

  double u_l2error = exact_sol.l2_error(system_name, "u");
  double u_h1error = exact_sol.h1_error(system_name, "u");

  double v_l2error = exact_sol.l2_error(system_name, "v");
  double v_h1error = exact_sol.h1_error(system_name, "v");

  double p_l2error = exact_sol.l2_error(system_name, "p");
  double p_h1error = exact_sol.h1_error(system_name, "p");

  double T_l2error = exact_sol.l2_error(system_name, "T");
  double T_h1error = exact_sol.h1_error(system_name, "T");

  double wN_l2error = exact_sol.l2_error(system_name, "w_N");
  double wN_h1error = exact_sol.h1_error(system_name, "w_N");

  double wN2_l2error = exact_sol.l2_error(system_name, "w_N2");
  double wN2_h1error = exact_sol.h1_error(system_name, "w_N2");

  double wCN_l2error = exact_sol.l2_error(system_name, "w_CN");
  double wCN_h1error = exact_sol.h1_error(system_name, "w_CN");
  
  double w_l2error = 0.0, 
         w_h1error = 0.0;

  if( (es->get_mesh()).mesh_dimension() == 3 )
    {
      w_l2error = exact_sol.l2_error(system_name, "w");
      w_h1error = exact_sol.h1_error(system_name, "w");
    }

  int return_flag = 0;

  // This is the tolerance of the iterative linear solver so
  // it's unreasonable to expect anything better than this.
  double tol = 5.0e-10;
  
  if( u_l2error > tol   || u_h1error > tol   ||
      v_l2error > tol   || v_h1error > tol   ||
      w_l2error > tol   || w_h1error > tol   ||
      p_l2error > tol   || p_h1error > tol   ||
      T_l2error > tol   || T_h1error > tol   ||
      wN_l2error > tol  || wN_h1error > tol  ||
      wN2_l2error > tol || wN2_h1error > tol ||
      wCN_l2error > tol || wCN_h1error > tol )
    {
      return_flag = 1;

      std::cout << "Tolerance exceeded for thermally driven flow test." << std::endl
		<< "tolerance     = " << tol << std::endl
		<< "u l2 error    = " << u_l2error << std::endl
		<< "u h1 error    = " << u_h1error << std::endl
		<< "v l2 error    = " << v_l2error << std::endl
		<< "v h1 error    = " << v_h1error << std::endl
		<< "w l2 error    = " << w_l2error << std::endl
		<< "w h1 error    = " << w_h1error << std::endl
		<< "p l2 error    = " << p_l2error << std::endl
		<< "p h1 error    = " << p_h1error << std::endl
		<< "T l2 error    = " << T_l2error << std::endl
		<< "T h1 error    = " << T_h1error << std::endl
		<< "w_N l2 error  = " << wN_l2error << std::endl
		<< "w_N h1 error  = " << wN_h1error << std::endl
		<< "w_N2 l2 error = " << wN2_l2error << std::endl
		<< "w_N2 h1 error = " << wN2_h1error << std::endl
                << "w_CN l2 error = " << wCN_l2error << std::endl
                << "w_CN h1 error = " << wCN_h1error << std::endl;
    }

  return return_flag;
}

Real initial_values( const Point& p, const Parameters &params, 
		     const std::string& , const std::string& unknown_name )
{
  Real value = 0.0;

  if( unknown_name == "w_N2" )
    value = params.get<Real>("w_N2");

  else if( unknown_name == "w_N" )
    value = params.get<Real>("w_N");

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
