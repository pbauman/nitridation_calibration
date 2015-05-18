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
#include <full_model_composition.h>
#include <full_model_likelihood.h>

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

    NitridationCalibration::FullModelComposition<QUESO::GslVector,QUESO::GslMatrix>
      full_model(argc,argv,*env,model_input);

    QUESO::GslVector guess(full_model.get_model().param_space().zeroVector());

    unsigned int guess_size = model_input.vector_variable_size("Optimizer/guess");

    if( guess_size != guess.sizeLocal() )
      {
        if( env->fullRank() == 0 )
          {
            std::cerr << "Error: Initial guess size mismatch!" << std::endl
                      << "input guess size = " << guess_size << std::endl
                      << "sip guess size   = " << guess.sizeLocal() << std::endl;
          }
        env.reset();
        MPI_Finalize();
        return 1;
      }

    for( unsigned int i = 0; i < guess_size; i++ )
      {
        guess[i] = model_input("Optimizer/guess", 0.0, i);
        if( env->fullRank() == 0 )
          {
            std::cout << "guess["<<i<<"] = " << guess[i] << std::endl;
          }
      }

    //************************************************
    // Solve SIP and get back posterior RV
    //************************************************
    //sip->solve();
  }

  MPI_Finalize();

#else

  std::cout << "Must have linked against a valid QUESO installation for this program to run." << std::endl;

#endif //NITCAL_HAVE_QUESO

  return 0;
}
