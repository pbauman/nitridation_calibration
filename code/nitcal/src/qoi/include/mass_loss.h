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

#ifndef NITCAL_MASS_LOSS_H
#define NITCAL_MASS_LOSS_H

// C++
#include <set>
#include <vector>

// GRINS
#include "grins/qoi_base.h"

// Antioch
#include "antioch/vector_utils.h"
#include "antioch/chemical_mixture.h"

namespace GRINS
{
  class ReactingLowMachNavierStokesBase;
}

namespace NitridationCalibration
{
  class MassLoss : public GRINS::QoIBase
  {
  public:

    MassLoss( const GetPot& input );
    ~MassLoss();

    virtual libMesh::AutoPtr<libMesh::DifferentiableQoI> clone();

    virtual void side_qoi( DiffContext& context, const QoISet& qoi_indices );

    virtual void init( const GetPot& input, const GRINS::MultiphysicsSystem& system );

  protected:

    //! List of boundary ids for which we want to compute this QoI
    std::set<libMesh::boundary_id_type> _bc_ids;

    std::vector<GRINS::VariableIndex> _species_vars;

    GRINS::VariableIndex _T_var;

    GRINS::ReactingLowMachNavierStokesBase* _physics;

    unsigned int _CN_index;

    //! Scales mass flux to mass loss in kg
    libMesh::Real _factor;

    Antioch::ChemicalMixture<libMesh::Real>* _chem_mixture;

  private:

    MassLoss();

  };

} // end namespace NitridationCalibration

#endif // NITCAL_MASS_LOSS_H
