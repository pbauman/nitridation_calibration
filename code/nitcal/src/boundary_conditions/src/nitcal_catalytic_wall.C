//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "nitcal_catalytic_wall.h"

namespace NitridationCalibration
{

  CatalyticWall::CatalyticWall( const GetPot& input, 
				GRINS::VariableIndex T_var, 
				const std::vector<GRINS::VariableIndex>& species_vars,
				GRINS::VariableIndex species, 
				unsigned int species_index,
				const GRINS::ChemicalMixture& chem_mixture,
				Real gamma )
    : GRINS::NeumannFuncObj(),
      _T_var(T_var),
      _species_vars(species_vars),
      _p0( input("Physics/"+GRINS::reacting_low_mach_navier_stokes+"/p0", 0.0 ) ),
      _s_var(species),
      _s_index(species_index),
      _chem_mixture(chem_mixture),
      _gamma(gamma)
  {
    return;
  }

  CatalyticWall::~CatalyticWall()
  {
    return;
  }

  Real CatalyticWall::normal_value( const libMesh::FEMContext& context, const unsigned int qp )
  {
    // Get temperature at the quadrature point
    const Real T = context.side_value(_T_var, qp );
    
    const unsigned int n_species = _chem_mixture.species_list().size();
    std::vector<Real> mass_fractions( n_species );
    for( unsigned int s = 0 ; s < n_species; s++ )
      {
	mass_fractions[s] = context.side_value(_species_vars[s], qp );
      }

    const Real rho = _p0/(_chem_mixture.R(mass_fractions)*T);

    return omega_dot( rho*mass_fractions[_s_var], _chem_mixture.R(_s_index),
		      T, _chem_mixture.M(_s_index) );
  }
  
  Real CatalyticWall::normal_derivative( const libMesh::FEMContext& /*context*/, 
					 const unsigned int /*qp*/ )
  {
    libmesh_not_implemented();
    return 0.0;
  }

  Real CatalyticWall::normal_derivative( const libMesh::FEMContext& /*context*/, 
					 const unsigned int /*qp*/, 
					 const GRINS::VariableIndex /*jac_var*/ )
  {
    libmesh_not_implemented();
    return 0.0;
  }

} //namespace NitridationCalibration
