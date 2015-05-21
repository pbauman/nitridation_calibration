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
#include <surrogate_model_composition.h>
#include <surrogate_model_likelihood.h>
#include <queso_sip_interface.h>

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
  std::string model_inputfile = argv[1];
  std::string QUESO_input = argv[2];

  //************************************************
  // Initialize environments
  //************************************************
  MPI_Init(&argc,&argv);

  {
    boost::scoped_ptr<QUESO::FullEnvironment>
      env( new QUESO::FullEnvironment( MPI_COMM_WORLD,
                                       QUESO_input.c_str(),
                                       "",
                                       NULL ) );


    GetPot model_input( model_inputfile );

    NitridationCalibration::SurrogateModelComposition<QUESO::GslVector,QUESO::GslMatrix>
      surrogate_model(*env,model_input);

    // Read in points at which to evaluate the model
    std::ifstream input;
    std::vector<double> gamma_N_data;
    gamma_N_data.reserve(195000);
    std::vector<double> gamma_CN_data;
    gamma_CN_data.reserve(195000);

    input.open("datapoints.dat");

    std::cout << "Input is " << input.good() << std::endl;
    std::cout << "Reading Data!" << std::endl;
    while( input.good() )
      {
        double gCN, gN;

        input >> gCN >> gN;

        gamma_CN_data.push_back(gCN);
        gamma_N_data.push_back(gN);
      }

    std::cout << "gamma N size = " << gamma_N_data.size() << std::endl;
    std::cout << "gamma CN size = " << gamma_CN_data.size() << std::endl;
    input.close();
    std::cout << "Done reading data!" << std::endl;

    QUESO::GslVector params(surrogate_model.get_model().param_space().zeroVector());

    QUESO::GslVector output(surrogate_model.get_observations());

    std::ofstream ofile;
    ofile.open("modelpoints.dat");

    for( unsigned int i = 0; i < gamma_N_data.size(); i++)
      {
        params[0] = gamma_CN_data[i];
        params[1] = gamma_N_data[i];

        surrogate_model.compute_values( params, output);

        for( unsigned int j = 0; j < output.sizeGlobal(); j++ )
          ofile << output[j] << " ";

        ofile << std::endl;
      }

    ofile.close();

    MPI_Finalize();
  }

#else

  std::cout << "Must have linked against a valid QUESO installation for this program to run." << std::endl;

#endif //NITCAL_HAVE_QUESO

  return 0;
}
