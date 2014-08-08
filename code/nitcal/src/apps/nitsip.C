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

#include "nitcal_config.h"

// C++
#include <string>
#include <iostream>

// Boost
#include "boost/scoped_ptr.hpp"

// NitCal
#include "constant_gamma_cn_sip.h"
#include "constant_gamma_n_constant_gamma_cn_sip.h"
#include "arrhenius_gamma_n_constant_gamma_cn_sip.h"
#include "power_gamma_n_constant_gamma_cn_sip.h"
#include "arrhenius_gamma_n_arrhenius_gamma_cn_sip.h"
#include "power_gamma_n_arrhenius_gamma_cn_sip.h"
#include "arrhenius_gamma_n_power_gamma_cn_sip.h"
#include "power_gamma_n_power_gamma_cn_sip.h"

#ifdef NITCAL_HAVE_QUESO
// QUESO
#include "queso/GslVector.h"
#include "queso/GslMatrix.h"
#endif // HAVE_QUESO

int main(int argc, char* argv[])
{
#ifdef NITCAL_HAVE_QUESO
  // Check command line count.
  if( argc < 3 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify SIP and QUESO input files." << std::endl;
      exit(1); // TODO: something more sophisticated for parallel runs?
    }

  // libMesh input file should be first argument
  std::string sip_input_filename = argv[1];
  std::string QUESO_input = argv[2];

  //************************************************
  // Initialize environments
  //************************************************
  MPI_Init(&argc,&argv);

  QUESO::FullEnvironment* env = new QUESO::FullEnvironment(MPI_COMM_WORLD,
							   QUESO_input.c_str(),
							   "",
							   NULL );

  //************************************************
  // Need artifical block here because libMesh needs 
  // to call PetscFinalize before we call 
  // MPI_Finalize.
  //************************************************
  {
    GetPot sip_input(sip_input_filename);
    std::string sip_type = sip_input( "InverseProblem/sip_type", "DIE!" );

    boost::scoped_ptr<NitridationCalibration::StatisticalInverseProblemBase<QUESO::GslVector,QUESO::GslMatrix> > sip(NULL);

    if( sip_type == std::string("constant_gamma_cn") )
      {
        sip.reset( new NitridationCalibration::ConstantGammaCNSIP<QUESO::GslVector,QUESO::GslMatrix>( env, "multilevel", argc, argv, sip_input_filename ) );
      }
    else if( sip_type == std::string("constant_gamma_n_constant_gamma_cn") )
      {
        sip.reset( new NitridationCalibration::ConstantGammaNConstantGammaCNSIP<QUESO::GslVector,QUESO::GslMatrix>( env, "multilevel", argc, argv, sip_input_filename ) );
      }
    else if( sip_type == std::string("arrhenius_gamma_n_constant_gamma_cn") )
      {
        sip.reset( new NitridationCalibration::ArrheniusGammaNConstantGammaCNSIP<QUESO::GslVector,QUESO::GslMatrix>( env, "multilevel", argc, argv, sip_input_filename ) );
      }
    else if( sip_type == std::string("power_gamma_n_constant_gamma_cn") )
      {
        sip.reset( new NitridationCalibration::PowerGammaNConstantGammaCNSIP<QUESO::GslVector,QUESO::GslMatrix>( env, "multilevel", argc, argv, sip_input_filename ) );
      }
    else if( sip_type == std::string("arrhenius_gamma_n_arrhenius_gamma_cn") )
      {
        sip.reset( new NitridationCalibration::ArrheniusGammaNArrheniusGammaCNSIP<QUESO::GslVector,QUESO::GslMatrix>( env, "multilevel", argc, argv, sip_input_filename ) );
      }
    else if( sip_type == std::string("power_gamma_n_arrhenius_gamma_cn") )
      {
        sip.reset( new NitridationCalibration::PowerGammaNArrheniusGammaCNSIP<QUESO::GslVector,QUESO::GslMatrix>( env, "multilevel", argc, argv, sip_input_filename ) );
      }
    else if( sip_type == std::string("arrhenius_gamma_n_power_gamma_cn") )
      {
        sip.reset( new NitridationCalibration::ArrheniusGammaNPowerGammaCNSIP<QUESO::GslVector,QUESO::GslMatrix>( env, "multilevel", argc, argv, sip_input_filename ) );
      }
    else if( sip_type == std::string("power_gamma_n_power_gamma_cn") )
      {
        sip.reset( new NitridationCalibration::PowerGammaNPowerGammaCNSIP<QUESO::GslVector,QUESO::GslMatrix>( env, "multilevel", argc, argv, sip_input_filename ) );
      }
    else
      {
        std::cerr << "Error: Invalid SIP type! Found " << sip_type << std::endl;
        delete env;
        MPI_Finalize();
        return 1;
      }

    //************************************************
    // Solve SIP and get back posterior RV
    //************************************************
    sip->solve();
  }

  //************************************************
  // Finalize environments
  //************************************************
  delete env;

  MPI_Finalize();

#else

  std::cout << "Must have linked against a valid QUESO installation for this program to run." << std::endl;

#endif //NITCAL_HAVE_QUESO
  
  return 0;
}
