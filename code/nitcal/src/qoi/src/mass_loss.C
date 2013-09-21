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

// This class
#include "mass_loss.h"

// GRINS
#include "grins/multiphysics_sys.h"
#include "grins/reacting_low_mach_navier_stokes_base.h"
#include "grins/variable_name_defaults.h"
#include "grins/math_constants.h"

// libMesh
#include "libmesh/quadrature.h"

namespace NitridationCalibration
{
  MassLoss::MassLoss( const GetPot& input )
    : QoIBase(),
      _physics(NULL),
      _chem_mixture(NULL)
  {
    this->assemble_qoi_sides = true;
    this->assemble_qoi_elements = false;

    unsigned int n_species = input.vector_variable_size("Physics/Chemistry/species");
    std::vector<std::string> species_list(n_species);

    for( unsigned int s = 0; s < n_species; s++ )
      {
        species_list[s] = input( "Physics/Chemistry/species", "DIE!", s );
      }

    _chem_mixture =  new Antioch::ChemicalMixture<libMesh::Real>( species_list );

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

    libMesh::Real radius = input("QoI/MassLoss/radius", 0.0 );

    libMesh::Real length = input("QoI/MassLoss/length", 0.0 );

    libMesh::Real delta_t = input("QoI/MassLoss/delta_t", 0.0 );

    // delta_t*total_area = delta_t*( 2*(cap areas) + length*surface area )
    this->_factor = delta_t*( GRINS::Constants::two_pi*radius*radius + GRINS::Constants::two_pi*radius*length );

    return;
  }

  MassLoss::~MassLoss()
  {
    return;
  }

  libMesh::AutoPtr<libMesh::DifferentiableQoI> MassLoss::clone()
  {
    return libMesh::AutoPtr<libMesh::DifferentiableQoI>( new MassLoss( *this ) );
  }

  void MassLoss::init( const GetPot& input, const GRINS::MultiphysicsSystem& system )
  {
    std::tr1::shared_ptr<GRINS::Physics> base_physics = system.get_physics( GRINS::reacting_low_mach_navier_stokes );
    _physics = libmesh_cast_ptr<GRINS::ReactingLowMachNavierStokesBase* >( base_physics.get() );
    
    // Grab temperature variable index
    std::string T_var_name = input("Physics/VariableNames/Temperature",
				   GRINS::T_var_name_default);

    this->_T_var = system.variable_number(T_var_name);

    const unsigned int n_species = _physics->n_species();

    _species_vars.resize( n_species );
    for( unsigned int s = 0; s < n_species; s++ )
      {
	std::string var_name = "w_"+_chem_mixture->chemical_species()[s]->species();
	_species_vars[s] = system.variable_number( var_name );
      }

    _CN_index = _chem_mixture->active_species_name_map().find(std::string("CN"))->second;
    
    return;
  }

  void MassLoss::side_qoi( libMesh::DiffContext& context, const libMesh::QoISet& )
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

	    unsigned int n_qpoints = c.get_side_qrule().n_points();

	    const std::vector<libMesh::Point>& normals = side_fe->get_normals();
	    
	    libMesh::Number& qoi = c.get_qois()[0];

            const unsigned int n_species = _physics->n_species();

            std::vector<libMesh::Real> Y, D;
            Y.resize(n_species);
            D.resize(n_species);

            for (unsigned int qp = 0; qp != n_qpoints; qp++)
              {
                const libMesh::Real T =  c.side_value(_T_var, qp);
                
                for( unsigned int s = 0; s < n_species; s++ )
                  {
                    c.side_value( _species_vars[s], qp, Y[s] );
                  }
                
                const libMesh::Real p0 = _physics->get_p0_steady_side(c, qp);
                
                const libMesh::Real R = _chem_mixture->R( Y );

                const libMesh::Real rho = _physics->rho(T, p0, R);

                const libMesh::Real cp = _physics->cp_mix( T, Y );

                const libMesh::Real k = _physics->k( T, Y );

                _physics->D( rho, cp, k, D );
		
		libMesh::Gradient grad_Y;
		c.side_gradient( _species_vars[_CN_index], qp, grad_Y );

		qoi += _factor*rho*D[_CN_index]*grad_Y*normals[qp]*JxW[qp];

	      } // quadrature loop

	  } // end check on boundary id
      }

    return;
  }

} // end namespace NitridationCalibration
