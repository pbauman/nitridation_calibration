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

// Antioch
#include "antioch/vector_utils.h"

// This class
#include "average_N_mole_fraction.h"

// GRINS
#include "grins/multiphysics_sys.h"
#include "grins/reacting_low_mach_navier_stokes_base.h"
#include "grins/variable_name_defaults.h"
#include "grins/math_constants.h"
#include "grins/assembly_context.h"

// libMesh
#include "libmesh/quadrature.h"

namespace NitridationCalibration
{
  AverageNMoleFraction::AverageNMoleFraction( const std::string& qoi_name )
    : QoIBase(qoi_name)
  {
    return;
  }

  AverageNMoleFraction::~AverageNMoleFraction()
  {
    return;
  }

  GRINS::QoIBase* AverageNMoleFraction::clone() const
  {
    return new AverageNMoleFraction( *this );
  }

  void AverageNMoleFraction::init( const GetPot& input,
                                   const GRINS::MultiphysicsSystem& system )
  {
    const unsigned int n_species = input.vector_variable_size("Physics/Chemistry/species");

    std::vector<std::string> species_list(n_species);

    for( unsigned int s = 0; s < n_species; s++ )
      {
        species_list[s] = input( "Physics/Chemistry/species", "DIE!", s );
      }

    _chem_mixture =  new Antioch::ChemicalMixture<libMesh::Real>( species_list );

    // Read boundary ids for which we want to compute
    int num_bcs =  input.vector_variable_size("QoI/AverageNMoleFraction/bc_ids");

    if( num_bcs <= 0 )
      {
	std::cerr << "Error: Must specify at least one boundary id to compute"
		  << " average N mole fraction." << std::endl
		  << "Found: " << num_bcs << std::endl;
	libmesh_error();
      }

    for( int i = 0; i < num_bcs; i++ )
      {
	_bc_ids.insert( input("QoI/AverageNMoleFraction/bc_ids", -1, i ) );
      }

    libMesh::Real radius = input("QoI/AverageNMoleFraction/channel_radius", 0.0 );

    // Factor is 2*pi (for integrating dtheta)/(pi*radius^2) = 2/r^2
    this->_factor = 2.0/(radius*radius);

    _species_vars.resize( n_species );
    for( unsigned int s = 0; s < n_species; s++ )
      {
	std::string var_name = "w_"+_chem_mixture->chemical_species()[s]->species();
	_species_vars[s] = system.variable_number( var_name );
      }

    _N_index = _chem_mixture->active_species_name_map().find(std::string("N"))->second;

    return;
  }

  void AverageNMoleFraction::init_context( GRINS::AssemblyContext& context )
  {
    libMesh::FEBase* N_fe;

    context.get_side_fe<libMesh::Real>(this->_species_vars[_N_index], N_fe);

    N_fe->get_phi();
    N_fe->get_JxW();

    return;
  }

  void AverageNMoleFraction::side_qoi( GRINS::AssemblyContext& context,
                                       const unsigned int qoi_index )
  {
    for( std::set<libMesh::boundary_id_type>::const_iterator id = _bc_ids.begin();
	 id != _bc_ids.end(); id++ )
      {
	if( context.has_side_boundary_id( (*id) ) )
	  {
	    FEBase* side_fe;
	    context.get_side_fe<libMesh::Real>(this->_species_vars[_N_index], side_fe);

	    const std::vector<libMesh::Real> &JxW = side_fe->get_JxW();

	    unsigned int n_qpoints = context.get_side_qrule().n_points();
	    
            const std::vector<libMesh::Point>& qpoint = side_fe->get_xyz();

	    libMesh::Number& qoi = context.get_qois()[qoi_index];

            const unsigned int n_species = _chem_mixture->n_species();

            std::vector<libMesh::Real> Y;
            Y.resize(n_species);

            for (unsigned int qp = 0; qp != n_qpoints; qp++)
              {
                const libMesh::Real r = qpoint[qp](0);

                for( unsigned int s = 0; s < n_species; s++ )
                  {
                    context.side_value( _species_vars[s], qp, Y[s] );
                  }
                
                const libMesh::Real M = _chem_mixture->M( Y );

                const libMesh::Real Y_N = Y[_N_index];

                libMesh::Real X_N = _chem_mixture->X<libMesh::Real>( _N_index, M, Y_N );

		qoi += _factor*X_N*r*JxW[qp];

	      } // quadrature loop

	  } // end check on boundary id
      }

    return;
  }

} // end namespace NitridationCalibration
