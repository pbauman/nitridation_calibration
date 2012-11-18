//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// NitCal - Nitridation Calibration 
//
// Copyright (C) 2010-2012 The PECOS Development Team
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
#include "getpot.h"
#include "libmesh_common.h"
#include "point.h"
#include "function_base.h"


namespace NitridationCalibration
{

  class TubeTempBC : public libMesh::FunctionBase<Real>
  {
  public:

    TubeTempBC( const GetPot& input );
    virtual ~TubeTempBC();

    virtual libMesh::AutoPtr<libMesh::FunctionBase<Real> > clone() const;

    virtual Real operator()(const libMesh::Point& p, const Real time=0.);

    virtual void operator()(const libMesh::Point& p, const Real time, 
			    libMesh::DenseVector<Real>& output);

  protected:

    Real linear_interp( const Real x ) const;
    Real spline_interp( const Real ) const
    { libmesh_not_implemented(); }

    std::vector<Real> _wall_tc_locs;
    std::vector<Real> _wall_temps;

  };

} // namespace NitridationCalibration

#endif //NITCAL_TUBE_TWALL_H
