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

// This class
#include "average_N_mole_fraction.h"

namespace NitridationCalibration
{
  AverageNMoleFraction::AverageNMoleFraction( const std::string& qoi_name )
    : QoIBase(qoi_name)
  {
    return;
  }

  AverageNMoleFraction::~AverageNMoleFraction()
  {
    return;
  }

  GRINS::QoIBase* AverageNMoleFraction::clone() const
  {
    return new AverageNMoleFraction( *this );
  }

  void AverageNMoleFraction::init( const GetPot& input,
                                   const GRINS::MultiphysicsSystem& system )
  {
    const unsigned int n_species = input.vector_variable_size("Physics/Chemistry/species");

    std::vector<std::string> species_list(n_species);

    for( unsigned int s = 0; s < n_species; s++ )
      {
        species_list[s] = input( "Physics/Chemistry/species", "DIE!", s );
      }

    // Read boundary ids for which we want to compute
    int num_bcs =  input.vector_variable_size("QoI/MassLoss/bc_ids");

    if( num_bcs <= 0 )
      {
	std::cerr << "Error: Must specify at least one boundary id to compute"
		  << " mass loss." << std::endl
		  << "Found: " << num_bcs << std::endl;
	libmesh_error();
      }

    for( int i = 0; i < num_bcs; i++ )
      {
	_bc_ids.insert( input("QoI/MassLoss/bc_ids", -1, i ) );
      }

    libMesh::Real radius = input("QoI/MassLoss/radius", 0.0 );

    _species_vars.resize( n_species );
    for( unsigned int s = 0; s < n_species; s++ )
      {
	std::string var_name = "w_"+_chem_mixture->chemical_species()[s]->species();
	_species_vars[s] = system.variable_number( var_name );
      }

    _N_index = _chem_mixture->active_species_name_map().find(std::string("N"))->second;

    return;
  }

} // end namespace NitridationCalibration
