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

#ifndef NITCAL_CATALYTIC_WALL_H
#define NITCAL_CATALYTIC_WALL_H

// libMesh
#include "getpot.h"

// GRINS
#include "grins/math_constants.h"
#include "grins/grins_physics_names.h"
#include "grins/neumann_func_obj.h"
#include "grins/chemical_mixture.h"

namespace NitridationCalibration
{

  class CatalyticWall : public GRINS::NeumannFuncObj
  {
  public:

    CatalyticWall( const GetPot& input, 
		   GRINS::VariableIndex T_var, 
		   const std::vector<GRINS::VariableIndex>& species_vars,
		   GRINS::VariableIndex species, 
		   unsigned int species_index,
		   const GRINS::ChemicalMixture& chem_mixture,
		   Real gamma );

    virtual ~CatalyticWall();

    //! Returns the value of the implemented Neumann boundary condition
    /*! This will leverage the FEMContext to get variable values and derivatives through the
      side_value, side_gradient, etc. interfaces, for each quadrature point qp. */
    virtual libMesh::Point value( const libMesh::FEMContext&, const unsigned int )
    { libmesh_error(); return libMesh::Point(); }

    //! Returns the value of the implemented Neumann boundary condition
    /*! This will leverage the FEMContext to get variable values and derivatives through the
      side_value, side_gradient, etc. interfaces, for each quadrature point qp. 
      Returns the normal component of the Neumann value. Only to be used when flux vector is
      formulated implicitly in terms of normal component. By default, does nothing since
      it's only applicable in special cases. */
    virtual Real normal_value( const libMesh::FEMContext& context, const unsigned int qp );

    //! Returns the derivative with respect to the primary variable of the implemented Neumann boundary condition.
    /*! This will leverage the FEMContext to get variable values and derivatives through the
      side_value, side_gradient, etc. interfaces, for each quadrature point qp. */
    virtual libMesh::Point derivative( const libMesh::FEMContext&, const unsigned int )
    { libmesh_error(); return libMesh::Point(); }

    //! If needed, returns the derivative with respect to other variables in the system.
    /*! By default, does nothing. User should reimplement is this is needed.
      This will leverage the FEMContext to get variable values and derivatives through the
      side_value, side_gradient, etc. interfaces, for each quadrature point qp. */
    virtual libMesh::Point derivative( const libMesh::FEMContext&, 
				       const unsigned int, 
				       const GRINS::VariableIndex )
    { libmesh_error(); return libMesh::Point(); }

    //! Returns the derivative with respect to the primary variable of the implemented Neumann boundary condition.
    /*! This will leverage the FEMContext to get variable values and derivatives through the
      side_value, side_gradient, etc. interfaces, for each quadrature point qp.
      Returns the normal component of the Neumann value. Only to be used when flux vector is
      formulated implicitly in terms of normal component. By default, does nothing since
      it's only applicable in special cases. */
    virtual Real normal_derivative( const libMesh::FEMContext& context, const unsigned qp );


    //! If needed, returns the derivative with respect to other variables in the system.
    /*! By default, does nothing. User should reimplement is this is needed.
      This will leverage the FEMContext to get variable values and derivatives through the
      side_value, side_gradient, etc. interfaces, for each quadrature point qp.
      Returns the normal component of the Neumann value. Only to be used when flux vector is
      formulated implicitly in terms of normal component. By default, does nothing since
      it's only applicable in special cases. */
    virtual Real normal_derivative( const libMesh::FEMContext& context, const unsigned int qp, 
				    const GRINS::VariableIndex jac_var );

  protected:

    inline
    Real omega_dot( Real rho_s, Real R_s, Real T, Real M_s ) const
    {
      return rho_s*_gamma*std::sqrt( (R_s*T)/(2*GRINS::Constants::pi*M_s) );
    }
   
    const GRINS::VariableIndex _T_var;
    std::vector<GRINS::VariableIndex> _species_vars;

    const Real _p0;

    const GRINS::VariableIndex _s_var;
    unsigned int _s_index;

    GRINS::ChemicalMixture _chem_mixture;

    const Real _gamma;

  private:

    CatalyticWall();

  };

}// namespace NitridationCalibration

#endif //NITCAL_CATALYTIC_WALL_H
