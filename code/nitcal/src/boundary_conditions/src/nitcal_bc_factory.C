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

#include "nitcal_bc_factory.h"

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
    
    GRINS::DBCContainer cont;
    cont.add_var_name( "T" );
    cont.add_bc_id( 4 );
    std::tr1::shared_ptr<libMesh::FunctionBase<libMesh::Real> > tube_temp( new TubeTempBC( _input ) );
    cont.set_func( tube_temp );


    std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > dbc;

    dbc.insert( std::make_pair(GRINS::reacting_low_mach_navier_stokes,cont) );

    return dbc;
  }

} // NitridationCalibration
