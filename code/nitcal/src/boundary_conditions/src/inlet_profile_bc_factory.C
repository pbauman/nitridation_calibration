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

// This class
#include "inlet_profile_bc_factory.h"

// NitCal
#include "inlet_profile.h"

// Antioch
#include "antioch/chemical_mixture.h"

// GRINS
#include "grins/multiphysics_sys.h"
#include "grins/species_mass_fracs_fe_variables.h"
#include "grins/variable_warehouse.h"
#include "grins/variables_parsing.h"
#include "grins/string_utils.h"

namespace NitridationCalibration
{
  libMesh::UniquePtr<libMesh::FunctionBase<libMesh::Real> >
  InletProfileBCFactory::build_func( const GetPot& input,
                                     GRINS::MultiphysicsSystem& /*system*/,
                                     std::vector<std::string>& var_names,
                                     const std::string& section )
  {
    if( section.find("Velocity") == std::string::npos )
      libmesh_error_msg("ERROR: TubeTempBC is only valid for Velocity vars!");

    if( !input.have_variable( section+"/mdot" ) )
      libmesh_error_msg("Error: Must specify input mass flow rate mdot!");

    // We're expecting this in mg/s
    libMesh::Real mdot = input( section+"/mdot", 0.0 );

    // convert to kg/s
    mdot /= 1.0e6;

    if( !input.have_variable( section+"/p_inlet" ) )
      libmesh_error_msg("Error: Must specify "+section+"/p_inlet");

    libMesh::Real P = input(section+"/p_inlet", 600);

    std::string T_input_str = this->T_wall_input(section);
    libMesh::Real T = input(T_input_str, 0);

    // This only makes sense for SpeciesMassFractionsFEVariables
    // in the VariableWarehouse. This call will error out if it's not there.
    const GRINS::SpeciesMassFractionsFEVariables& species_fe_var =
      GRINS::GRINSPrivate::VariableWarehouse::get_variable_subclass<GRINS::SpeciesMassFractionsFEVariables>
      (GRINS::VariablesParsing::species_mass_fractions_section());

    unsigned int n_species = var_names.size();
    std::vector<std::string> mole_fracs_names(n_species);
    std::vector<std::string> species_names(n_species);

    const std::string& prefix = species_fe_var.prefix();
    for( unsigned int v = 0; v < n_species; v++ )
      {
        this->extract_species_name(var_names[v],prefix,species_names[v]);
        mole_fracs_names[v] = "X_"+species_names[v];
      }

    Antioch::ChemicalMixture<libMesh::Real> chem_mixture( species_names );

    std::vector<libMesh::Real> X( n_species );

    for( unsigned int s = 0; s < n_species; s++ )
      {
        std::string mf_input_str = this->mole_fracs_input(section,mole_fracs_names,s);
        X[s] = input(mf_input_str, 0.0 );
      }

    libMesh::Real M = 0.0;
    for( unsigned int s = 0; s < n_species; s++ )
      M += X[s]*chem_mixture.M(s);

    libMesh::Real R = Antioch::Constants::R_universal<libMesh::Real>()/M;

    libMesh::Real rho = P/(R*T);

    if( !input.have_variable( section+"/r0" ) )
      libmesh_error_msg("Error: Must specify input radius r0!");

    libMesh::Real r0 = input( section+"/r0", 0.0 );

    libMesh::UniquePtr<libMesh::FunctionBase<libMesh::Real> >
      func( new InletProfile(r0,mdot,rho) );

    return func;
  }

  std::string InletProfileBCFactory::T_wall_input(const std::string& section) const
  {
    std::vector<std::string> tokens;
    GRINS::StringUtilities::split_string(section,"/",tokens);

    return tokens[0]+"/"+tokens[1]+"/Temperature/T";
  }

  std::string InletProfileBCFactory::mole_fracs_input(const std::string& section,
                                                      std::vector<std::string>& mole_fracs_names,
                                                      unsigned int s) const
  {
    std::vector<std::string> tokens;
    GRINS::StringUtilities::split_string(section,"/",tokens);

    return tokens[0]+"/"+tokens[1]+"/SpeciesMassFractions/"+mole_fracs_names[s];
  }

  void InletProfileBCFactory::extract_species_name( const std::string& var_name,
                                                    const std::string& prefix,
                                                    std::string& species_name ) const
  {
    std::vector<std::string> split_name;
    GRINS::StringUtilities::split_string(var_name,prefix,split_name);

    // split_string won't add the prefix, since it was used as the delimiter, so
    // split_name should just have the lingering species name
    libmesh_assert_equal_to(split_name.size(), 1);

    species_name = split_name[0];
  }

} // end namespace NitridationCalibration
