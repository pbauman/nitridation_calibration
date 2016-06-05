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

#include "nitcal_config.h"

#include <iostream>

// GRINS
#include "grins/simulation.h"
#include "grins/simulation_builder.h"

// libMesh
#include "libmesh/parallel.h"

// NitCal
#include "qoi_factory.h"
#include "tube_twall.h"

// Function for getting initial temperature field
libMesh::Real initial_values( const libMesh::Point& p, const libMesh::Parameters &params,
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

  GetPot command_line(argc,argv);

  // Initialize libMesh library.
  libMesh::LibMeshInit libmesh_init(argc, argv);

  GRINS::SimulationBuilder sim_builder;

  GRINS::SharedPtr<GRINS::QoIFactory> qoi_factory( new NitridationCalibration::QoIFactory );

  sim_builder.attach_qoi_factory( qoi_factory );

  GRINS::Simulation grins( libMesh_inputfile,
                           command_line,
			   sim_builder,
                           libmesh_init.comm());

  //FIXME: We need to move this to within the Simulation object somehow...
  std::string restart_file = libMesh_inputfile( "restart-options/restart_file", "none" );



  if( restart_file == "none" )
    {
      // Asssign initial temperature value
      std::string system_name = libMesh_inputfile( "screen-options/system_name", "GRINS" );
      GRINS::SharedPtr<libMesh::EquationSystems> es = grins.get_equation_system();
      const libMesh::System& system = es->get_system(system_name);

      libMesh::Parameters &params = es->parameters;

      libMesh::Real& w_N2 = params.set<libMesh::Real>( "w_N2" );
      w_N2 = libMesh_inputfile( "BoundaryConditions/Inlet/SpeciesMassFractions/X_N2", 0.0);

      libMesh::Real& w_N = params.set<libMesh::Real>( "w_N" );
      w_N = libMesh_inputfile( "BoundaryConditions/Inlet/SpeciesMassFractions/X_N", 0.0);


      GRINS::SharedPtr<NitridationCalibration::TubeTempBC> wall_temp;

      std::string tc_loc_str("BoundaryConditions/OuterWall/Temperature/tc_locs");
      unsigned int tc_size =
      libMesh_inputfile.vector_variable_size(tc_loc_str);

      std::vector<libMesh::Real> wall_tc_locs(tc_size);

      for( unsigned int i = 0; i < tc_size; i++ )
        wall_tc_locs[i] = libMesh_inputfile(tc_loc_str, 0.0, i );

      if( !libMesh_inputfile.have_variable("BoundaryConditions/OuterWall/Temperature/wall_temps") )
        libmesh_error_msg("ERROR: Could not find BoundaryConditions/OuterWall/Temperature/wall_temps in input!");

      std::string wall_temp_str("BoundaryConditions/OuterWall/Temperature/wall_temps");
      unsigned int temp_size = libMesh_inputfile.vector_variable_size(wall_temp_str);

      if( temp_size != tc_size )
        libmesh_error_msg("Error: Must be same number of wall temp locations and wall temps.");

      std::vector<libMesh::Real> wall_temps(temp_size);

      for( unsigned int i = 0; i < temp_size; i++ )
        wall_temps[i] = libMesh_inputfile(wall_temp_str, 0.0, i );

      wall_temp.reset( new NitridationCalibration::TubeTempBC(wall_tc_locs,wall_temps) );
      GRINS::SharedPtr<NitridationCalibration::TubeTempBC>& dummy = params.set<GRINS::SharedPtr<NitridationCalibration::TubeTempBC> >( "wall_temp" );
      dummy = wall_temp;

      system.project_solution( initial_values, NULL, params );
    }

  grins.run();

  return 0;
}

libMesh::Real initial_values( const libMesh::Point& p, const libMesh::Parameters &params,
                              const std::string& , const std::string& unknown_name )
{
  libMesh::Real value = 0.0;

  if( unknown_name == "w_N2" )
    value = params.get<libMesh::Real>("w_N2");

  else if( unknown_name == "w_N" )
    value = params.get<libMesh::Real>("w_N");

  else if( unknown_name == "T" )
    {
      value = (*params.get<GRINS::SharedPtr<NitridationCalibration::TubeTempBC> >( "wall_temp" ))(p);
      //value = 1200.0;
    }

  else if( unknown_name == "u" )
    //value = 40.0;
  value = 0.0;

  else
    value = 0.0;

  return value;
}
