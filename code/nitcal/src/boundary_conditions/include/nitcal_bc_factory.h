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

#ifndef NITCAL_BC_FACTORY_H
#define NITCAL_BC_FACTORY_H

#include <string>

// libMesh
#include "getpot.h"
#include "equation_systems.h"
#include "system.h"
#include "zero_function.h"

// GRINS
#include "bc_types.h"
#include "bc_factory.h"

// NitridationCalibration
#include "nitcal_catalytic_wall.h"
#include "tube_twall.h"

namespace NitridationCalibration
{
  //! Object for constructing NitridationCalibration specific boundary condition function objects.
  class BoundaryConditionsFactory : public GRINS::BoundaryConditionsFactory
  {
  public:
    
    BoundaryConditionsFactory( const GetPot& input );

    virtual ~BoundaryConditionsFactory();
    
    //! Builds the NitridationCalibration::
    virtual std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > build_dirichlet( );

    //! Builds the NitridationCalibration::
    virtual std::map< GRINS::PhysicsName, GRINS::NBCContainer > build_neumann( libMesh::EquationSystems& es );

  protected:

    const GetPot& _input;

  }; // class BoundaryConditionsFactory

} // namespace NitridationCalibration

#endif //NITCAL_BC_FACTORY_H
