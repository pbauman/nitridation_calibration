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

// This class
#include "mass_loss.h"

// GRINS
#include "grins/multiphysics_sys.h"
#include "grins/ideal_gas_mixture.h"
#include "grins/reacting_low_mach_navier_stokes.h"
#include "grins/grins_kinetics.h"
#include "grins/cea_thermo.h"
#include "grins/constant_transport.h"
#include "grins/cached_values.h"
#include "grins/variable_name_defaults.h"

// libMesh
#include "libmesh/quadrature.h"

namespace NitridationCalibration
{

  template<class Mixture>
  MassLoss<Mixture>::MassLoss( const GetPot& input )
    : QoIBase()
  {
    this->assemble_qoi_sides = true;
    this->assemble_qoi_elements = false;

    // Read boundary ids for which we want to compute
    int num_bcs =  input.vector_variable_size("QoI/MassLoss/bc_ids");

    if( num_bcs <= 0 )
      {
	std::cerr << "Error: Must specify at least one boundary id to compute"
		  << " mass loss." << std::endl
		  << "Found: " << num_bcs << std::endl;
	libmesh_error();
      }

    for( int i = 0; i < num_bcs; i++ )
      {
	_bc_ids.insert( input("QoI/MassLoss/bc_ids", -1, i ) );
      }

    return;
  }

  template<class Mixture>
  MassLoss<Mixture>::~MassLoss()
  {
    return;
  }

  template<class Mixture>
  libMesh::AutoPtr<libMesh::DifferentiableQoI> MassLoss<Mixture>::clone()
  {
    return libMesh::AutoPtr<libMesh::DifferentiableQoI>( new MassLoss( *this ) );
  }

  template<class Mixture>
  void MassLoss<Mixture>::init( const GetPot& input, const GRINS::MultiphysicsSystem& system )
  {
    std::tr1::shared_ptr<GRINS::Physics> base_physics = system.get_physics( GRINS::reacting_low_mach_navier_stokes );
    _physics = libmesh_cast_ptr<const GRINS::ReactingLowMachNavierStokes<Mixture>* >( base_physics.get() );
    
    // Grab temperature variable index
    std::string T_var_name = input("Physics/VariableNames/Temperature",
				   GRINS::T_var_name_default);

    this->_T_var = system.variable_number(T_var_name);

    for( unsigned int s = 0; s < _physics->gas_mixture().chem_mixture().n_species(); s++ )
      {
	std::string var_name = "w_"+_physics->gas_mixture().chem_mixture().species_inverse_name_map().find(_physics->gas_mixture().chem_mixture().species_list()[s])->second;
	_species_vars[s] = system.variable_number( var_name );
      }

    _CN_index = _physics->gas_mixture().chem_mixture().active_species_name_map().find(std::string("CN"))->second;
    
    return;
  }

  template<class Mixture>
  void MassLoss<Mixture>::side_qoi( libMesh::DiffContext& context, const libMesh::QoISet& )
  {
    libMesh::FEMContext &c = libmesh_cast_ref<libMesh::FEMContext&>(context);

    for( std::set<libMesh::boundary_id_type>::const_iterator id = _bc_ids.begin();
	 id != _bc_ids.end(); id++ )
      {
	if( c.has_side_boundary_id( (*id) ) )
	  {
	    FEBase* side_fe;
	    c.get_side_fe<libMesh::Real>(this->_species_vars[_CN_index], side_fe);

	    const std::vector<libMesh::Real> &JxW = side_fe->get_JxW();

	    unsigned int n_qpoints = (c.get_side_qrule())->n_points();

	    const std::vector<libMesh::Point>& normals = side_fe->get_normals();
	    
	    libMesh::Number& qoi = c.elem_qoi[0];

	    GRINS::CachedValues cache;
	    
	    // Build up cache
	    {
	      std::vector<libMesh::Real> T;
	      T.resize(n_qpoints);

	      std::vector<std::vector<libMesh::Real> > Y;
	      Y.resize( n_qpoints );

	      std::vector<libMesh::Real> rho;
	      rho.resize(n_qpoints);

	      // Build up cache
	      for (unsigned int qp = 0; qp != n_qpoints; qp++)
		{
		  T[qp] =  c.side_value(_T_var, qp);

		  for( unsigned int s = 0; s < _physics->gas_mixture().chem_mixture().n_species(); s++ )
		    {
		      Y[qp].resize(_physics->gas_mixture().chem_mixture().n_species());

		      c.side_value( _species_vars[s], qp, Y[qp][s] );
		    }

		  const libMesh::Real p0 = _physics->get_p0_steady_side(c, qp);

		  rho[qp] = _physics->rho(T[qp], p0, Y[qp]);
		}

	      cache.set_values( GRINS::Cache::TEMPERATURE, T );
	      cache.set_vector_values( GRINS::Cache::MASS_FRACTIONS, Y );
	      cache.set_values( GRINS::Cache::MIXTURE_DENSITY, rho );

	      std::vector<libMesh::Real> cp;
	      cp.resize(n_qpoints);

	      for (unsigned int qp = 0; qp != n_qpoints; qp++)
		{
		  cp[qp] = _physics->gas_mixture().cp( cache, qp );
		}
	      cache.set_values( GRINS::Cache::MIXTURE_SPECIFIC_HEAT_P, cp );
	      

	      std::vector<std::vector<libMesh::Real> > D;
	      D.resize(n_qpoints);
	      for (unsigned int qp = 0; qp != n_qpoints; qp++)
		{
		  for( unsigned int s = 0; s < _physics->gas_mixture().chem_mixture().n_species(); s++ )
		    {
		      D[qp].resize(_physics->gas_mixture().chem_mixture().n_species());
		      _physics->gas_mixture().D( cache, qp, D[qp] );
		    }
		}

	      cache.set_vector_values( GRINS::Cache::DIFFUSION_COEFFS, D );
	    }

	    // Loop over quadrature points  
	    for (unsigned int qp = 0; qp != n_qpoints; qp++)
	      {
		const libMesh::Real& rho = cache.get_cached_values(GRINS::Cache::MIXTURE_DENSITY)[qp];
		const libMesh::Real& D_CN = cache.get_cached_vector_values(GRINS::Cache::DIFFUSION_COEFFS)[qp][_CN_index];
		
		libMesh::Gradient grad_Y_CN;
		c.side_gradient( _species_vars[_CN_index], qp, grad_Y_CN );

		qoi += rho*D_CN*grad_Y_CN*normals[qp]*JxW[qp];

	      } // quadrature loop

	  } // end check on boundary id
      }

    return;
  }

  //Instantiate
  template class MassLoss< GRINS::IdealGasMixture<GRINS::CEAThermodynamics,GRINS::ConstantTransport,GRINS::Kinetics> >;

} // end namespace NitridationCalibration
