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

#ifndef NITCAL_INLET_PROFILE_H
#define NITCAL_INLET_PROFILE_H

// C++
#include <vector>

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/point.h"
#include "libmesh/function_base.h"


namespace NitridationCalibration
{

  class InletProfile : public libMesh::FunctionBase<libMesh::Real>
  {
  public:

    InletProfile( libMesh::Real r0,
                  libMesh::Real mdot,
                  libMesh::Real rho );

    virtual ~InletProfile(){};

    virtual libMesh::UniquePtr<libMesh::FunctionBase<libMesh::Real> > clone() const;

    virtual libMesh::Real operator()( const libMesh::Point& p, const libMesh::Real time = 0. );

    virtual libMesh::Real operator()( const libMesh::Point& p, const libMesh::Real time = 0. ) const;

    virtual void operator()( const libMesh::Point& p, const libMesh::Real time,
                             libMesh::DenseVector<libMesh::Real>& output );

  protected:

    libMesh::Real parabolic_profile( const libMesh::Real r ) const;

    const libMesh::Real _r0_sq;

    libMesh::Real _C;

  private:

    InletProfile();

  };

} // namespace NitridationCalibration

#endif //NITCAL_TUBE_TWALL_H
