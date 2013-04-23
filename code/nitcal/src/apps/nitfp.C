//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// NitCal - Nitridation Calibration 
//
// Copyright (C) 2012-2013 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
#include "nitcal_config.h"

#include <iostream>

// GRINS
#include "grins/simulation.h"
#include "grins/simulation_builder.h"

// libMesh
#include "libmesh/parallel.h"

// NitCal
#include "nitcal_bc_factory.h"
#include "qoi_factory.h"

// Function for getting initial temperature field
Real initial_values( const Point& p, const Parameters &params, 
		     const std::string& system_name, const std::string& unknown_name );

int main(int argc, char* argv[])
{
  // Check command line count.
  if( argc < 2 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify libMesh input file." << std::endl;
      exit(1); // TODO: something more sophisticated for parallel runs?
    }

  // libMesh input file should be first argument
  std::string libMesh_input_filename = argv[1];
  
  // Create our GetPot object.
  GetPot libMesh_inputfile( libMesh_input_filename );

  // Initialize libMesh library.
  LibMeshInit libmesh_init(argc, argv);
 
  GRINS::SimulationBuilder sim_builder;

  std::tr1::shared_ptr<GRINS::BoundaryConditionsFactory> bc_factory( new NitridationCalibration::BoundaryConditionsFactory(libMesh_inputfile) );

  sim_builder.attach_bc_factory( bc_factory );

  std::tr1::shared_ptr<GRINS::QoIFactory> qoi_factory( new NitridationCalibration::QoIFactory );

  sim_builder.attach_qoi_factory( qoi_factory );
  
  GRINS::Simulation grins( libMesh_inputfile,
			   sim_builder );

  //FIXME: We need to move this to within the Simulation object somehow...
  std::string restart_file = libMesh_inputfile( "restart-options/restart_file", "none" );

  std::tr1::shared_ptr<NitridationCalibration::TubeTempBC> wall_temp;
  
  if( restart_file == "none" )
    {
      // Asssign initial temperature value
      std::string system_name = libMesh_inputfile( "screen-options/system_name", "GRINS" );
      std::tr1::shared_ptr<libMesh::EquationSystems> es = grins.get_equation_system();
      const libMesh::System& system = es->get_system(system_name);
      
      Parameters &params = es->parameters;

      Real& w_N2 = params.set<Real>( "w_N2" );
      w_N2 = libMesh_inputfile( "Physics/ReactingLowMachNavierStokes/bound_species_1", 0.0, 0 );
      
      Real& w_N = params.set<Real>( "w_N" );
      w_N = libMesh_inputfile( "Physics/ReactingLowMachNavierStokes/bound_species_1", 0.0, 1 );

      wall_temp.reset( new NitridationCalibration::TubeTempBC( libMesh_inputfile ) );
      std::tr1::shared_ptr<NitridationCalibration::TubeTempBC>& dummy = params.set<std::tr1::shared_ptr<NitridationCalibration::TubeTempBC> >( "wall_temp" );
      dummy = wall_temp;

      system.project_solution( initial_values, NULL, params );
    }

  grins.run();

  return 0;
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
