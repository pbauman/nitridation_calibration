//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// NitCal - Nitridation Calibration 
//
// Copyright (C) 2012 The PECOS Development Team
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

#include "nitcal_bc_factory.h"

namespace NitridationCalibration
{
  BoundaryConditionsFactory::BoundaryConditionsFactory( const GetPot& input )
    : GRINS::BoundaryConditionsFactory(), 
      _input(input),
      _T_var_name(input( "Physics/VariableNames/Temperature", "T" )),
      _system_name( input( "screen-options/system_name", "NitCal" ) )
  {
    unsigned int n_species = input.vector_variable_size("Physics/Chemistry/species");
    _species_names.resize( n_species, "DIE!" );

    for( unsigned int i = 0; i < n_species; i++ )
      {
        _species_names[i] = input( "Physics/Chemistry/species", "DIE!", i );
      }

    
    _chem_mixture = new GRINS::ChemicalMixture( _species_names );

    return;
  }

  BoundaryConditionsFactory::~BoundaryConditionsFactory()
  {
    delete _chem_mixture;
    return;
  }

  std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > BoundaryConditionsFactory::build_dirichlet( )
  {
    std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > dbc;

    return dbc;
  }

  std::map< GRINS::PhysicsName, GRINS::NBCContainer > BoundaryConditionsFactory::build_neumann( libMesh::EquationSystems& es )
  {
    const libMesh::System& system = es.get_system(_system_name);

    const GRINS::VariableIndex T_var = system.variable_number( _T_var_name );
    const GRINS::VariableIndex N_var = system.variable_number( "w_N" );
    const GRINS::VariableIndex CN_var = system.variable_number( "w_CN" );

    std::vector<GRINS::VariableIndex> species_vars( _species_names.size() );
    for( unsigned int s = 0; s < _species_names.size(); s++ )
      {
	species_vars[s] = system.variable_number( "w_"+_species_names[s] );
      }

    // Hard coded from input file for the moment.
    const unsigned int N_index = 1;
    const unsigned int CN_index = 3;

    const Real gamma_CN = _input( "Physics/BoundaryConditions/CatalyticWall/gamma_CN", -1.0 );
    const Real gamma_N = -gamma_CN;

    std::tr1::shared_ptr<GRINS::NeumannFuncObj> N_wall( new CatalyticWall(_input, T_var, species_vars, N_var, N_index, *_chem_mixture, gamma_N) );

    std::tr1::shared_ptr<GRINS::NeumannFuncObj> CN_wall( new CatalyticWall(_input,T_var, species_vars, CN_var, CN_index, *_chem_mixture, gamma_CN) );

    GRINS::NeumannBCsMap N_map;
    GRINS::NeumannBCsMap CN_map;

    N_map.insert( GRINS::NBCMapPair(N_var, N_wall) );
    CN_map.insert( GRINS::NBCMapPair(CN_var, CN_wall) );

    const GRINS::BoundaryID wall_id = _input( "Physics/BoundaryConditions/CatalyticWall/catalytic_wall_id", -1 );

    GRINS::NBCContainer nbc_container;
    nbc_container.insert( GRINS::NBCContainerPair( wall_id, N_map ) );
    nbc_container.insert( GRINS::NBCContainerPair( wall_id, CN_map ) );

    std::map< GRINS::PhysicsName, GRINS::NBCContainer > nbc;

    nbc.insert( std::make_pair(GRINS::reacting_low_mach_navier_stokes, nbc_container) );

    return nbc;
  }

} // NitridationCalibration
