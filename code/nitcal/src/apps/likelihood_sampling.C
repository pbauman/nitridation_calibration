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
#include <iostream>
#include <iomanip>
#include <fstream>

// NitCal
#include "nitcal_config.h"
#include "gamma_cn_sip.h"
#include "gamma_cn_likelihood.h"

#ifdef NITCAL_HAVE_QUESO
// QUESO
#include "uqGslVector.h"
#include "uqGslMatrix.h"
#endif // HAVE_QUESO

int main(int argc, char* argv[])
{
#ifdef NITCAL_HAVE_QUESO
  // Check command line count.
  if( argc < 3 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify libMesh and QUESO input files." << std::endl;
      exit(1); // TODO: something more sophisticated for parallel runs?
    }

  // libMesh input file should be first argument
  std::string libMesh_input_filename = argv[1];
  std::string QUESO_input = argv[2];

  //************************************************
  // Initialize environments
  //************************************************
  MPI_Init(&argc,&argv);

  uqFullEnvironmentClass* env = new uqFullEnvironmentClass(MPI_COMM_WORLD,
							   QUESO_input.c_str(),
							   "",
							   NULL );
  {
    NitridationCalibration::GammaCNSIP<uqGslVectorClass,uqGslMatrixClass> sip( env, "multilevel",
									       argc, argv,
									       libMesh_input_filename );

    const uqBaseScalarFunctionClass<uqGslVectorClass,uqGslMatrixClass>& raw_likelihood = sip.get_likelihood_func();

    const NitridationCalibration::GammaCNLikelihood<uqGslVectorClass,uqGslMatrixClass>& likelihood = 
      libmesh_cast_ref<const NitridationCalibration::GammaCNLikelihood<uqGslVectorClass,uqGslMatrixClass>& >( raw_likelihood );

    std::ofstream values;
    std::ofstream errors;

    std::vector<double> param( 1, 0.01);

    for( unsigned int s = 0; s < 200; s++ )
      {
        double likelihood_value = 0.0;
        likelihood.update_parameters(param);

        try
          {
            likelihood_value = likelihood.evaluate_likelihood();
          }
        catch(...)
          {
            if( env->fullRank() == 0 )
              {
                std::cout << "Caught exception for evaluating likelihood for parameter = " << param[0] << std::endl;
                errors.open( "error_params.dat", std::ios::app );
                errors << std::scientific << std::setprecision(16) << param[0] << std::endl;
                errors.close();
              }

            param[0] += 0.01;
            continue;
          }

        if( env->fullRank() == 0 )
          {
            values.open( "likelihood_values.dat", std::ios::app );
            values << std::scientific << std::setprecision(16)
                   << param[0] << " " << likelihood_value << std::endl;
            values.close();
          }
        
        param[0] += 0.01;
      }
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
