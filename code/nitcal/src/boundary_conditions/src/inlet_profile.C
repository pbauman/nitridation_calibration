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

// This class
#include "inlet_profile.h"

// GRINS
#include "grins/math_constants.h"

// Antioch
#include "antioch/physical_constants.h"
#include "antioch/chemical_mixture.h"

namespace NitridationCalibration
{
  InletProfile::InletProfile( libMesh::Real r0,
                              libMesh::Real mdot,
                              libMesh::Real rho )
    : _r0_sq(r0*r0),
      _C(mdot/rho*2.0/_r0_sq/GRINS::Constants::pi)
  {}

  libMesh::UniquePtr<libMesh::FunctionBase<libMesh::Real> > InletProfile::clone() const
  {
    return libMesh::UniquePtr<libMesh::FunctionBase<libMesh::Real> >( new InletProfile(*this) );
  }

  libMesh::Real InletProfile::operator()(const libMesh::Point& p, const libMesh::Real /*time*/ )
  {
    const libMesh::Real& r = p(0);

    return this->parabolic_profile( r );
  }

  libMesh::Real InletProfile::operator()(const libMesh::Point& p, const libMesh::Real /*time*/ ) const
  {
    const libMesh::Real& r = p(0);

    return this->parabolic_profile( r );
  }

  void InletProfile::operator()( const libMesh::Point& p, const libMesh::Real time,
                                 libMesh::DenseVector<libMesh::Real>& output )
  {
    for( unsigned int i = 0; i < output.size(); i++ )
      {
	output(i) = (*this)(p,time);
      }

    return;
  }

  libMesh::Real InletProfile::parabolic_profile( const libMesh::Real r ) const
  {
    return _C*(1.0 - (r*r)/_r0_sq);
  }

} // end namespace NitridationCalibration
