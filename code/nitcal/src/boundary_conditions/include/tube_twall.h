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

#ifndef NITCAL_TUBE_TWALL_H
#define NITCAL_TUBE_TWALL_H

// C++
#include <vector>

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/point.h"
#include "libmesh/function_base.h"


namespace NitridationCalibration
{

  class TubeTempBC : public libMesh::FunctionBase<libMesh::Real>
  {
  public:

    TubeTempBC( const GetPot& input );
    virtual ~TubeTempBC();

    virtual libMesh::AutoPtr<libMesh::FunctionBase<libMesh::Real> > clone() const;

    virtual libMesh::Real operator()(const libMesh::Point& p, const libMesh::Real time=0.);

    virtual libMesh::Real operator()(const libMesh::Point& p, const libMesh::Real time=0.) const;

    virtual void operator()(const libMesh::Point& p, const libMesh::Real time, 
			    libMesh::DenseVector<libMesh::Real>& output);

  protected:

    libMesh::Real linear_interp( const libMesh::Real x ) const;
    libMesh::Real spline_interp( const libMesh::Real ) const
    { libmesh_not_implemented(); }

    std::vector<libMesh::Real> _wall_tc_locs;
    std::vector<libMesh::Real> _wall_temps;

  private:

    TubeTempBC();

  };

} // namespace NitridationCalibration

#endif //NITCAL_TUBE_TWALL_H
