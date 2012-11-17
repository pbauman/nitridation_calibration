//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// VAR Simulator - 2-D Axisymmetric Finite Element Formulation 
//
// Copyright (C) 2011 The VAR Development Team
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
// $Id: var_bc_factory.C 94 2012-04-26 22:09:58Z pbauman $
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "var_bc_factory.h"

VAR::BoundaryConditionsFactory::BoundaryConditionsFactory( const GetPot& input )
  : GRINS::BoundaryConditionsFactory( input ),
    _input( input ),
    _T_var_name( input( "Physics/VariableNames/Temperature", "T" ) ),
    _u_r_var_name( input("Physics/VariableNames/r_velocity", "u_r" ) ),
    _u_z_var_name( input("Physics/VariableNames/z_velocity", "u_z" ) ),
    _system_name( input( "screen-options/system_name", "VAR" ) ),
    _top_bcid( input("mesh-options/top_bcid", 2) ), /*Default is set for libMesh built 2D grid*/
    _side_bcid( input("mesh-options/side_bcid", 1) ), /*Default is set for libMesh built 2D grid*/
    _bottom_bcid( input("mesh-options/bottom_bcid", 0) ), /*Default is set for libMesh built 2D grid*/
    _inflow( std::tr1::shared_ptr<libMesh::FunctionBase<Number> >() )
{
  return;
}

VAR::BoundaryConditionsFactory::~BoundaryConditionsFactory()
{
  return;
}

std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > VAR::BoundaryConditionsFactory::build_dirichlet()
{
  /* We can't guarantee the build order, so we check if the pointer 
     has been set or not and build it if it hasn't be set. */
  if( !_inflow )
    _inflow.reset( new VAR::InflowTop( _input ) );

  if( !_side_bottom_flow )
    _side_bottom_flow.reset( new VAR::SideBottomFlow( _input ) );

  GRINS::DBCContainer cont;
  cont.add_var_name( "z_vel" );
  cont.add_bc_id( _top_bcid );
  cont.set_func( _inflow );

  GRINS::DBCContainer cont2;
  cont2.add_var_name( "z_vel" );
  cont2.add_bc_id( _side_bcid );
  cont2.add_bc_id( _bottom_bcid );
  cont2.set_func( _side_bottom_flow );

  GRINS::DBCContainer cont3;
  cont3.add_var_name( "r_vel" );
  cont3.add_bc_id( _side_bcid );
  cont3.add_bc_id( _bottom_bcid );
  std::tr1::shared_ptr<libMesh::FunctionBase<Number> > zero( new  ZeroFunction<Number> );
  cont3.set_func( zero );
  
  std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > mymap;

  mymap.insert( std::pair<GRINS::PhysicsName, GRINS::DBCContainer >(GRINS::axisymmetric_incomp_navier_stokes,  cont) );

  mymap.insert( std::pair<GRINS::PhysicsName, GRINS::DBCContainer >(GRINS::axisymmetric_incomp_navier_stokes,  cont2) );

  mymap.insert( std::pair<GRINS::PhysicsName, GRINS::DBCContainer >(GRINS::axisymmetric_incomp_navier_stokes,  cont3) );

  return mymap;
}

std::map< std::string, GRINS::NBCContainer > VAR::BoundaryConditionsFactory::build_neumann( libMesh::EquationSystems& es )
{
  const libMesh::System& system = es.get_system(_system_name);

  const GRINS::VariableIndex T_var = system.variable_number( _T_var_name );

  /* We can't guarantee the build order, so we check if the pointer has been set or not
     and build it if it hasn't be set. */
  if( !_inflow )
    {
      _inflow.reset( new VAR::InflowTop( _input ) );
    }

  if( !_side_bottom_flow )
    _side_bottom_flow.reset( new VAR::SideBottomFlow( _input ) );

  std::tr1::shared_ptr<GRINS::NeumannFuncObj> hf_top( new VAR::HeatFluxTop( _input, T_var, _inflow ) );
  std::tr1::shared_ptr<GRINS::NeumannFuncObj> hf_side( new VAR::HeatFluxSide( _input, T_var ) );
  std::tr1::shared_ptr<GRINS::NeumannFuncObj> hf_bottom( new VAR::HeatFluxBottom( _input, T_var, _side_bottom_flow ) );

  GRINS::NeumannBCsMap nbc_map_top, nbc_map_side, nbc_map_bottom;
  nbc_map_top.insert( GRINS::NBCMapPair( T_var, hf_top ) );
  nbc_map_side.insert( GRINS::NBCMapPair( T_var, hf_side ) );
  nbc_map_bottom.insert( GRINS::NBCMapPair( T_var, hf_bottom ) );

  GRINS::NBCContainer nbc_container;
  nbc_container.insert( GRINS::NBCContainerPair( _top_bcid, nbc_map_top ) );
  nbc_container.insert( GRINS::NBCContainerPair( _side_bcid, nbc_map_side ) );
  nbc_container.insert( GRINS::NBCContainerPair( _bottom_bcid, nbc_map_bottom ) );

  std::map< std::string, GRINS::NBCContainer > nbcs;
  nbcs.insert( std::pair< std::string, GRINS::NBCContainer >( "AxisymmetricHeatTransfer", nbc_container ) );

  return nbcs;
}
