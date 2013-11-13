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
#include "nitcal_bc_factory.h"

// NitCal
#include "inlet_profile.h"

namespace NitridationCalibration
{
  BoundaryConditionsFactory::BoundaryConditionsFactory( const GetPot& input )
    : GRINS::BoundaryConditionsFactory(), 
      _input(input),
      _T_var_name(input( "Physics/VariableNames/Temperature", "T" )),
      _system_name( input( "screen-options/system_name", "NitCal" ) )
  {
    return;
  }

  BoundaryConditionsFactory::~BoundaryConditionsFactory()
  {
    //delete _chem_mixture;
    return;
  }

  std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > BoundaryConditionsFactory::build_dirichlet( )
  {
    
    GRINS::DBCContainer T_cont;
    T_cont.add_var_name( "T" );
    T_cont.add_bc_id( 2 );
    std::tr1::shared_ptr<libMesh::FunctionBase<libMesh::Real> > tube_temp( new TubeTempBC( _input ) );
    T_cont.set_func( tube_temp );

    GRINS::DBCContainer v_cont;
    v_cont.add_var_name( "v" );
    v_cont.add_bc_id( 1 );
    std::tr1::shared_ptr<libMesh::FunctionBase<libMesh::Real> > v_inlet( new InletProfile( _input ) );
    v_cont.set_func( v_inlet );

    std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > dbc;

    dbc.insert( std::make_pair(GRINS::reacting_low_mach_navier_stokes, T_cont) );
    dbc.insert( std::make_pair(GRINS::reacting_low_mach_navier_stokes, v_cont) );

    return dbc;
  }

} // NitridationCalibration
