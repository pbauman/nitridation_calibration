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

#ifndef NITCAL_BC_FACTORY_H
#define NITCAL_BC_FACTORY_H

#include <string>

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/equation_systems.h"
#include "libmesh/system.h"
#include "libmesh/zero_function.h"

// GRINS
#include "grins/bc_types.h"
#include "grins/bc_factory.h"

// NitridationCalibration
//#include "nitcal_catalytic_wall.h"
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

  protected:

    const GetPot& _input;

    std::vector<std::string> _species_names;

    const std::string _T_var_name;
    const std::string _system_name;

  }; // class BoundaryConditionsFactory

} // namespace NitridationCalibration

#endif //NITCAL_BC_FACTORY_H
