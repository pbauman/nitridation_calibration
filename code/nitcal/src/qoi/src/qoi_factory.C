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
#include "qoi_factory.h"

// NitCal
#include "mass_loss.h"

// GRINS
#include "grins/ideal_gas_mixture.h"
#include "grins/grins_kinetics.h"
#include "grins/cea_thermo.h"
#include "grins/constant_transport.h"

namespace NitridationCalibration
{
  QoIFactory::QoIFactory()
    : GRINS::QoIFactory()
  {
    return;
  }

  QoIFactory::~QoIFactory()
  {
    return;
  }

  void QoIFactory::add_qoi( const GetPot& input,
			    const std::string& qoi_name,
			    std::tr1::shared_ptr<GRINS::QoIBase>& qoi )
  {
    if( qoi_name == "MassLoss" )
      {
	std::string chem_lib = input( "Physics/"+GRINS::reacting_low_mach_navier_stokes+"/chemistry_library", "cantera" );
	std::string thermo_lib = input( "Physics/"+GRINS::reacting_low_mach_navier_stokes+"/thermodynamics_library", "cantera" );
	std::string transport_lib = input( "Physics/"+GRINS::reacting_low_mach_navier_stokes+"/transport_library", "cantera" );

	if( chem_lib == "grins" && thermo_lib == "grins_cea" && transport_lib == "grins_constant" )
	  {
	    qoi.reset( new MassLoss<GRINS::IdealGasMixture< GRINS::CEAThermodynamics,GRINS::ConstantTransport,GRINS::Kinetics > >( input ) );
	  }
	else
	  {
	    std::cerr << "Invalid options for chemistry, thermo, and/or transport for MassLoss QoI:" << std::endl
		      << "requested chemistry lib = " << chem_lib << std::endl
		      << "requested thermo lib    = " << thermo_lib << std::endl
		      << "requested transport lib = " << transport_lib << std::endl;
	    libmesh_error();
	  }
      }
    
    else
      {
	GRINS::QoIFactory::add_qoi( input, qoi_name, qoi );
      }

    return;
  }

} // end namespace NitridationCalibration
