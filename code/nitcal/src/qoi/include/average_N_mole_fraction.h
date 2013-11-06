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

#ifndef NITCAL_AVERAGE_N_MOLE_FRACTION_H
#define NITCAL_AVERAGE_N_MOLE_FRACTION_H

// GRINS
#include "grins/qoi_base.h"

// Antioch
#include "antioch/vector_utils.h"
#include "antioch/chemical_mixture.h"

namespace NitridationCalibration
{

  class AverageNMoleFraction : public GRINS::QoIBase
  {
  public:

    AverageNMoleFraction( const std::string& qoi_name );

    virtual ~AverageNMoleFraction();

    virtual GRINS::QoIBase* clone() const;

    virtual bool assemble_on_interior() const;

    virtual bool assemble_on_sides() const;

    virtual void init( const GetPot& input, const GRINS::MultiphysicsSystem& system );

    virtual void init_context( GRINS::AssemblyContext& context );

    virtual void side_qoi( GRINS::AssemblyContext& context,
                           const unsigned int qoi_index );

  protected:

    std::set<libMesh::boundary_id_type> _bc_ids;

    std::vector<GRINS::VariableIndex> _species_vars;

    unsigned int _N_index;

    Antioch::ChemicalMixture<libMesh::Real>* _chem_mixture;

    libMesh::Real _factor;

  };

  inline
  bool AverageNMoleFraction::assemble_on_interior() const
  {
    return false;
  }

  inline
  bool AverageNMoleFraction::assemble_on_sides() const
  {
    return true;
  }

} // end namespace NitridationCalibration

#endif // NITCAL_AVERAGE_N_MOLE_FRACTION_H
