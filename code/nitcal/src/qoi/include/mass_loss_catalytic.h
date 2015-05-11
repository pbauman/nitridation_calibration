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
// $Id: mass_loss.h 42475 2013-11-06 18:59:50Z pbauman $
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef NITCAL_MASS_LOSS_CATALYTIC_H
#define NITCAL_MASS_LOSS_CATALYTIC_H

// C++
#include <set>
#include <vector>

// Antioch
#include "antioch/vector_utils_decl.h"

// GRINS
#include "grins/qoi_base.h"
#include "grins/antioch_chemistry.h"
#include "grins/gas_solid_catalytic_wall.h"

namespace GRINS
{
  class ReactingLowMachNavierStokesAbstract;
  class AssemblyContext;
}

namespace NitridationCalibration
{
  class MassLossCatalytic : public GRINS::QoIBase
  {
  public:

    MassLossCatalytic( const std::string& qoi_name );
    virtual ~MassLossCatalytic();

    virtual GRINS::QoIBase* clone() const;

    virtual bool assemble_on_interior() const;

    virtual bool assemble_on_sides() const;

    virtual void init( const GetPot& input, const GRINS::MultiphysicsSystem& system );

    virtual void init_context( GRINS::AssemblyContext& context );

    virtual void side_qoi( GRINS::AssemblyContext& context,
                           const unsigned int qoi_index );

    virtual void side_qoi_derivative( GRINS::AssemblyContext& context,
                                      const unsigned int qoi_index );

  protected:

    //! List of boundary ids for which we want to compute this QoI
    std::set<libMesh::boundary_id_type> _bc_ids;

    std::vector<GRINS::VariableIndex> _species_vars;

    GRINS::VariableIndex _T_var;

    GRINS::ReactingLowMachNavierStokesAbstract* _physics;

    unsigned int _N_index;
    unsigned int _CN_index;

    //! Scales mass flux to mass loss in kg
    libMesh::Real _factor;

    GRINS::AntiochChemistry* _chem_mixture;

    GRINS::GasSolidCatalyticWall<GRINS::AntiochChemistry>* _omega_dot;

  private:

    MassLossCatalytic();

  };
  
  inline
  bool MassLossCatalytic::assemble_on_interior() const
  {
    return false;
  }

  inline
  bool MassLossCatalytic::assemble_on_sides() const
  {
    return true;
  }

} // end namespace NitridationCalibration

#endif // NITCAL_MASS_LOSS_CATALYTIC_H
