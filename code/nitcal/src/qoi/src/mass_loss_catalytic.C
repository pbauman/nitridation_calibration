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
// $Id: mass_loss.C 42475 2013-11-06 18:59:50Z pbauman $
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// This class
#include "mass_loss_catalytic.h"

// GRINS
#include "grins/multiphysics_sys.h"
#include "grins/reacting_low_mach_navier_stokes_base.h"
#include "grins/variable_name_defaults.h"
#include "grins/math_constants.h"
#include "grins/assembly_context.h"
#include "grins/bc_handling_base.h"
#include "grins/reacting_low_mach_navier_stokes_bc_handling.h"

// libMesh
#include "libmesh/quadrature.h"

namespace NitridationCalibration
{
  MassLossCatalytic::MassLossCatalytic( const std::string& qoi_name )
    : QoIBase(qoi_name),
      _physics(NULL),
      _chem_mixture(NULL),
      _omega_dot(NULL)
  {
    return;
  }

  MassLossCatalytic::~MassLossCatalytic()
  {
    return;
  }

  GRINS::QoIBase* MassLossCatalytic::clone() const
  {
    return new MassLossCatalytic( *this );
  }

  void MassLossCatalytic::init( const GetPot& input, const GRINS::MultiphysicsSystem& system )
  {
    
    const unsigned int n_species = input.vector_variable_size("Physics/Chemistry/species");

    _chem_mixture =  new GRINS::AntiochChemistry(input);

    // Read boundary ids for which we want to compute
    int num_bcs =  input.vector_variable_size("QoI/MassLossCatalytic/bc_ids");

    if( num_bcs <= 0 )
      {
        std::cerr << "Error: Must specify at least one boundary id to compute"
                  << " mass loss." << std::endl
                  << "Found: " << num_bcs << std::endl;
        libmesh_error();
      }
    
     /*! \todo Need to generalize to more than 1 bc */
    if( num_bcs > 1 )
      {
        std::cerr << "Error: "+_qoi_name+" only support 1 bc_id at this time."
                  << std::endl
                  << "Found: " << num_bcs << std::endl;
        libmesh_error();
      }

    for( int i = 0; i < num_bcs; i++ )
      {
        _bc_ids.insert( input("QoI/"+_qoi_name+"/bc_ids", -1, i ) );
      }

    libMesh::Real delta_t = input("QoI/"+_qoi_name+"/delta_t", 0.0 );
    
    // convert from mins to seconds
    delta_t *= 60;

    // delta_t*\int d\theta
    this->_factor = delta_t*GRINS::Constants::two_pi;

    std::tr1::shared_ptr<GRINS::Physics> base_physics = system.get_physics( GRINS::reacting_low_mach_navier_stokes );
    _physics = libMesh::libmesh_cast_ptr<GRINS::ReactingLowMachNavierStokesBase* >( base_physics.get() );
    
    // Grab temperature variable index
    std::string T_var_name = input("Physics/VariableNames/Temperature",
                                   GRINS::T_var_name_default);

    this->_T_var = system.variable_number(T_var_name);

    _species_vars.resize( n_species );
    for( unsigned int s = 0; s < n_species; s++ )
      {
        std::string var_name = "w_"+_chem_mixture->species_name(s);
        _species_vars[s] = system.variable_number( var_name );
      }

    _N_index = _chem_mixture->species_index(std::string("N"));
    _CN_index = _chem_mixture->species_index(std::string("CN"));

    GRINS::BCHandlingBase* bc_handler_base = base_physics->get_bc_handler();

    GRINS::ReactingLowMachNavierStokesBCHandling<GRINS::AntiochChemistry>* bc_handler =
      libMesh::libmesh_cast_ptr<GRINS::ReactingLowMachNavierStokesBCHandling<GRINS::AntiochChemistry>*>(bc_handler_base);
    
    /*! \todo Need to generalize to more than 1 bc */
    GRINS::CatalyticWallBase<GRINS::AntiochChemistry>* base_wall_ptr  = bc_handler->get_catalytic_wall( (*_bc_ids.begin()) );
    _omega_dot = libMesh::libmesh_cast_ptr<GRINS::GasSolidCatalyticWall<GRINS::AntiochChemistry>*>( base_wall_ptr );

    return;
  }

  void MassLossCatalytic::init_context( GRINS::AssemblyContext& context )
  {
    libMesh::FEBase* CN_fe;

    context.get_side_fe<libMesh::Real>(this->_species_vars[_CN_index], CN_fe);

    CN_fe->get_phi();
    CN_fe->get_JxW();

    return;
  }

  void MassLossCatalytic::side_qoi( GRINS::AssemblyContext& context,
                                    const unsigned int qoi_index )
  {
    for( std::set<libMesh::boundary_id_type>::const_iterator id = _bc_ids.begin();
         id != _bc_ids.end(); id++ )
      {
        if( context.has_side_boundary_id( (*id) ) )
          {
            GRINS::VariableIndex CN_var = this->_species_vars[_CN_index];

            libMesh::FEBase* side_fe;
            context.get_side_fe<libMesh::Real>(CN_var, side_fe);

            const std::vector<libMesh::Real> &JxW = side_fe->get_JxW();

            unsigned int n_qpoints = context.get_side_qrule().n_points();

            const std::vector<libMesh::Point>& qpoint = side_fe->get_xyz();
            
            libMesh::Number& qoi = context.get_qois()[qoi_index];

            const unsigned int n_species = _physics->n_species();

            std::vector<libMesh::Real> Y;
            Y.resize(n_species);

            for (unsigned int qp = 0; qp != n_qpoints; qp++)
              {
                const libMesh::Real T =  context.side_value(_T_var, qp);
                
                const libMesh::Real r = qpoint[qp](0);

                for( unsigned int s = 0; s < n_species; s++ )
                  {
                    context.side_value( _species_vars[s], qp, Y[s] );
                  }
                
                const libMesh::Real p0 = _physics->get_p0_steady_side(context, qp);
                
                const libMesh::Real R = _chem_mixture->R_mix( Y );

                const libMesh::Real rho = _physics->rho(T, p0, R);

                const libMesh::Real mdot_C = _omega_dot->compute_reactant_solid_mass_consumption(rho,Y[_N_index],T);

                /*
                const libMesh::Real n_N = _chem_mixture->molar_density( _N_index, rho, Y[_N_index] )*1000;
                const libMesh::Real omega_dot_CN = _omega_dot->omega_dot(rho_N,T);
                std::cout << "rho_N = " << rho_N << std::endl;
                std::cout << "n_N = " << n_N << std::endl;
                std::cout << "T = " << T << std::endl;
                std::cout << "R_N = " << _chem_mixture->R(_N_index) << std::endl;
                std::cout << "v_N = " << std::sqrt( _chem_mixture->R(_N_index)*T/GRINS::Constants::two_pi) << std::endl;

                std::cout << "omega_dot = " << omega_dot_CN << std::endl;
                */

                qoi += _factor*mdot_C*r*JxW[qp];

              } // quadrature loop

          } // end check on boundary id
      }

    return;
  }

void MassLossCatalytic::side_qoi_derivative( GRINS::AssemblyContext& context,
                                             const unsigned int qoi_index )
  {
    for( std::set<libMesh::boundary_id_type>::const_iterator id = _bc_ids.begin();
         id != _bc_ids.end(); id++ )
      {
        if( context.has_side_boundary_id( (*id) ) )
          {
            GRINS::VariableIndex N_var = this->_species_vars[_N_index];

            const unsigned int n_T_dofs = context.get_dof_indices(_T_var).size();
            const unsigned int n_s_dofs = context.get_dof_indices(N_var).size();

            libMesh::FEBase* s_side_fe;
            context.get_side_fe<libMesh::Real>(N_var, s_side_fe);

            libMesh::FEBase* T_side_fe;
            context.get_side_fe<libMesh::Real>(_T_var, T_side_fe);

            const std::vector<std::vector<libMesh::Real> >& s_phi = s_side_fe->get_phi();

            const std::vector<std::vector<libMesh::Real> >& T_phi = T_side_fe->get_phi();

            const std::vector<libMesh::Real> &JxW = s_side_fe->get_JxW();

            unsigned int n_qpoints = context.get_side_qrule().n_points();

            const std::vector<libMesh::Point>& qpoint = s_side_fe->get_xyz();

            libMesh::DenseSubVector<libMesh::Number>& dQ_dT = context.get_qoi_derivatives(qoi_index, _T_var);

            libMesh::DenseSubVector<libMesh::Number>& dQ_dYN = context.get_qoi_derivatives(qoi_index, N_var);

            const unsigned int n_species = _physics->n_species();

            std::vector<libMesh::Real> Y;
            Y.resize(n_species);

            for (unsigned int qp = 0; qp != n_qpoints; qp++)
              {
                const libMesh::Real T =  context.side_value(_T_var, qp);
                
                const libMesh::Real r = qpoint[qp](0);

                for( unsigned int s = 0; s < n_species; s++ )
                  {
                    context.side_value( _species_vars[s], qp, Y[s] );
                  }
                
                const libMesh::Real p0 = _physics->get_p0_steady_side(context, qp);
                
                const libMesh::Real R = _chem_mixture->R_mix( Y );

                const libMesh::Real rho = _physics->rho(T, p0, R);

                const libMesh::Real dmdot_dT = _omega_dot->compute_reactant_solid_mass_consumption_dT(rho,Y[_N_index],T);

                const libMesh::Real dmdot_dYN = _omega_dot->compute_reactant_solid_mass_consumption_dYs(rho,Y,T);

                for( unsigned int i = 0; i != n_T_dofs; i++ )
                  {
                    dQ_dT(i) += _factor*dmdot_dT*T_phi[i][qp]*r*JxW[qp];
                  }

                for( unsigned int i = 0; i != n_s_dofs; i++ )
                  {
                    dQ_dYN(i) += _factor*dmdot_dYN*s_phi[i][qp]*r*JxW[qp];
                  }

              } // quadrature loop

          } // end check on boundary id
      }

    return;
  }

} // end namespace NitridationCalibration
