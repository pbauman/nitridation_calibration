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
    : _input(input)
  {
    return;
  }

  BoundaryConditionsFactory::~BoundaryConditionsFactory()
  {
    return;
  }

  std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > BoundaryConditionsFactory::build_dirichlet( )
  {
    std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > dbc;

    return dbc;
  }

  std::map< GRINS::PhysicsName, GRINS::NBCContainer > BoundaryConditionsFactory::build_neumann( libMesh::EquationSystems& es )
  {
    std::map< GRINS::PhysicsName, GRINS::NBCContainer > nbc;

    return nbc;
  }

} // NitridationCalibration
