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
#include "inlet_profile.h"

// GRINS
#include "grins/chemical_mixture.h"

namespace NitridationCalibration
{
  InletProfile::InletProfile( const GetPot& input )
    : _r0_sq( input( "BoundaryConditions/InletProfile/r0", 0.011 )*input( "BoundaryConditions/InletProfile/r0", 0.011 ) ),
      _C(-1.0)
  {
    if( !input.have_variable( "BoundaryConditions/InletProfile/Q" ) )
      {
        std::cerr << "Error: Must specify input volume flow rate Q!" << std::endl;
        libmesh_error();
      }

    // We're expecting this in sccm (cm^3/min)
    libMesh::Real Q = input( "BoundaryConditions/InletProfile/Q", 0.0 );

    // Convert Q to m^3/s
    Q /= (100*100*100*60);

    /* Now, covert to mass flux. According to 
       http://www.massflow-online.com/faqs/what-do-lnmin-lsmin-slm-and-sccm-stand-for
       conditions for sccm are 0 degrees C and 1 atm. */

    libMesh::Real P = 101325;
    libMesh::Real T = 273.15;

    unsigned int n_species = input.vector_variable_size( "Physics/Chemistry/species" );
    std::vector<std::string> species_names( n_species );
    for( unsigned int s = 0; s < n_species; s++ )
      {
	species_names[s] = input( "Physics/Chemistry/species", "DIE!", s );
      }

    GRINS::ChemicalMixture chem_mixture( species_names );

    /*
    std::vector<libMesh::Real> X( n_species );
    for( unsigned int s = 0; s < n_species; s++ )
      {
        X[s] = input( "Physics/ReactingLowMachNavierStokes/bound_species_1", 0.0 );
      }
    */

    //libMesh::Real R = chem_mixture.R_from_X( X ); ;
    //libMesh::Real R = chem_mixture.R( chem_mixture.species_name_map().find(std::string("N2"))->second );
    libMesh::Real R = chem_mixture.R( 0 );
    libMesh::Real rho = P/(R*T);

    Q *= rho;

    _C = 2.0*Q/_r0_sq/GRINS::Constants::pi;

    return;
  }
  
  InletProfile::~InletProfile()
  {
    return;
  }

  libMesh::AutoPtr<libMesh::FunctionBase<libMesh::Real> > InletProfile::clone() const
  {
    return libMesh::AutoPtr<libMesh::FunctionBase<libMesh::Real> >( new InletProfile(*this) );
  }

  libMesh::Real InletProfile::operator()(const libMesh::Point& p, const libMesh::Real /*time*/ )
  {
    const libMesh::Real& r = p(0);

    return this->parabolic_profile( r );
  }

  libMesh::Real InletProfile::operator()(const libMesh::Point& p, const libMesh::Real /*time*/ ) const
  {
    const libMesh::Real& r = p(0);

    return this->parabolic_profile( r );
  }

  void InletProfile::operator()( const libMesh::Point& p, const libMesh::Real time, 
                                 libMesh::DenseVector<libMesh::Real>& output )
  {
    for( unsigned int i = 0; i < output.size(); i++ )
      {
	output(i) = (*this)(p,time);
      }

    return;
  }

  libMesh::Real InletProfile::parabolic_profile( const libMesh::Real r ) const
  {
    return _C*(1.0 - (r*r)/_r0_sq);
  }

} // end namespace NitridationCalibration
