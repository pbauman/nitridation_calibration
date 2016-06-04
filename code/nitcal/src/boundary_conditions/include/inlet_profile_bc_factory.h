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

#ifndef NITCAL_INLET_PROFILE_BC_FACTORY_H
#define NITCAL_INLET_PROFILE_BC_FACTORY_H

// GRINS
#include "grins/dirichlet_bc_factory_function_base.h"

namespace GRINS
{
  class MultiphysicsSystem;
}

namespace NitridationCalibration
{
  class InletProfileBCFactory : public GRINS::DirichletBCFactoryFunctionBase<libMesh::FunctionBase<libMesh::Real> >
  {
  public:

    InletProfileBCFactory(const std::string& bc_type_name)
      : GRINS::DirichletBCFactoryFunctionBase<libMesh::FunctionBase<libMesh::Real> >(bc_type_name)
    {}

    virtual ~InletProfileBCFactory(){}

  protected:

    virtual libMesh::UniquePtr<libMesh::FunctionBase<libMesh::Real> >
    build_func( const GetPot& input,
                GRINS::MultiphysicsSystem& system,
                std::vector<std::string>& var_names,
                const std::string &section );

    std::string T_wall_input(const std::string& section) const;

    std::string mole_fracs_input(const std::string& section,
                                 std::vector<std::string>& mole_fracs_names,
                                 unsigned int s) const;

    void extract_species_name( const std::string& var_name,
                               const std::string& prefix,
                               std::string& species_name ) const;
  };

} // end namespace NitridationCalibration

#endif // NITCAL_INLET_PROFILE_BC_FACTORY_H
