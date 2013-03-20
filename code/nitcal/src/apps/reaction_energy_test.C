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

// C++
#include <string>
#include <vector>
#include <iostream>

// libMesh
#include "libmesh/getpot.h"

// GRINS
#include "grins/cea_thermo.h"

int main(int argc, char* argv[])
{

  std::vector<std::string> species_str_list;
  const unsigned int n_species = 4;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "N2" );
  species_str_list.push_back( "N" );
  species_str_list.push_back( "C" );
  species_str_list.push_back( "CN" );

  GRINS::ChemicalMixture chem_mixture( species_str_list );

  GRINS::CEAThermodynamics thermo( GetPot(), chem_mixture );

  const double P = 530.62;

  std::vector<double> T(3);
  T[0] = 1273.0;
  T[1] = 1000.0;
  T[2] = 600.0;

  for( std::vector<double>::const_iterator it = T.begin(); it != T.end(); ++it )
    {
      std::cout << "Delta_G(" << *it << ") = " << thermo.h_RT_minus_s_R(*it,2) +  thermo.h_RT_minus_s_R(*it,1) - thermo.h_RT_minus_s_R(*it,3) << std::endl;
    }

  return 0;
}
