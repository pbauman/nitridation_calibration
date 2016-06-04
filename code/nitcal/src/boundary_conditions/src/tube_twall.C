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
// General Public License for details more.
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

#include "tube_twall.h"

namespace NitridationCalibration
{

  TubeTempBC::TubeTempBC( const std::vector<libMesh::Real>& wall_tc_locs,
                          const std::vector<libMesh::Real>& wall_temps )
    : libMesh::FunctionBase<libMesh::Real>(),
    _wall_tc_locs(wall_tc_locs),
    _wall_temps(wall_temps)
  {}

  libMesh::UniquePtr<libMesh::FunctionBase<libMesh::Real> > TubeTempBC::clone() const
  {
    return libMesh::UniquePtr<libMesh::FunctionBase<libMesh::Real> >( new TubeTempBC(*this) );
  }

  libMesh::Real TubeTempBC::operator()( const libMesh::Point& p, const libMesh::Real )
  {
    const libMesh::Real& x = p(1);

    return this->linear_interp( x );
  }

  libMesh::Real TubeTempBC::operator()( const libMesh::Point& p, const libMesh::Real ) const
  {
    const libMesh::Real& x = p(1);

    return this->linear_interp( x );
  }

  void TubeTempBC::operator()( const libMesh::Point& p, const libMesh::Real time,
			       libMesh::DenseVector<libMesh::Real>& output )
  {
    for( unsigned int i = 0; i < output.size(); i++ )
      output(i) = (*this)(p,time);
  }

  libMesh::Real TubeTempBC::linear_interp( const libMesh::Real x ) const
  {
    // Find the bin
    unsigned int index = -1;

    // This is a stupid linear search. Should do a binary search.
    for( unsigned int i = 1; i < _wall_tc_locs.size(); i++ )
      {
	if( (x <= _wall_tc_locs[i]) ||
            (std::fabs( x -_wall_tc_locs[i])/_wall_tc_locs[i] < 1.0e-6)  )
	  {
	    index = i;
	    break;
	  }
      }

    libmesh_assert( index != static_cast<unsigned int>(-1) );

    // Interpolate
    const libMesh::Real x0 = _wall_tc_locs[index-1];
    const libMesh::Real x1 = _wall_tc_locs[index];

    const libMesh::Real T0 = _wall_temps[index-1];
    const libMesh::Real T1 = _wall_temps[index];

    return T0 + (T1-T0)/(x1-x0)*(x-x0);
  }

} // namespace NitridationCalibration
