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

// This class
#include "tube_wall_bc_factory.h"

// NitCal
#include "tube_twall.h"

// GRINS
#include "grins/multiphysics_sys.h"

namespace NitridationCalibration
{
  libMesh::UniquePtr<libMesh::FunctionBase<libMesh::Real> >
  TubeWallBCFactory::build_func( const GetPot& input,
                                 GRINS::MultiphysicsSystem& /*system*/,
                                 std::vector<std::string>& /*var_names*/,
                                 const std::string& section )
  {
    if( section.find("Temperature") == std::string::npos )
      libmesh_error_msg("ERROR: TubeTempBC is only valid for Temperature vars!");

    if( !input.have_variable(section+"/tc_locs") )
      libmesh_error_msg("ERROR: Could not find "+section+"/tc_locs in input!");

    std::string tc_loc_str(section+"/tc_locs");
    unsigned int tc_size =
      input.vector_variable_size(tc_loc_str);

    std::vector<libMesh::Real> wall_tc_locs(tc_size);

    for( unsigned int i = 0; i < tc_size; i++ )
      wall_tc_locs[i] = input(tc_loc_str, 0.0, i );

    if( !input.have_variable(section+"/wall_temps") )
      libmesh_error_msg("ERROR: Could not find "+section+"/wall_temps in input!");

    std::string wall_temp_str(section+"/wall_temps");
    unsigned int temp_size = input.vector_variable_size(wall_temp_str);

    if( temp_size != tc_size )
      libmesh_error_msg("Error: Must be same number of wall temp locations and wall temps.");

    std::vector<libMesh::Real> wall_temps(temp_size);

    for( unsigned int i = 0; i < temp_size; i++ )
      wall_temps[i] = input(wall_temp_str, 0.0, i );

    libMesh::UniquePtr<libMesh::FunctionBase<libMesh::Real> >
      func( new TubeTempBC(wall_tc_locs,wall_temps) );

    return func;
  }

  TubeWallBCFactory nitcal_tube_temp_bc_factory("nitridation_temp_bc");

} // end namespace NitridationCalibration
