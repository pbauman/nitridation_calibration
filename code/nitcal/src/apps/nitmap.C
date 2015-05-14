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
#include "queso/GslOptimizer.h"
#include "queso/OptimizerMonitor.h"
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

    // Initalize parameter vector
    QUESO::GslVector params( full_model.get_model().param_space().zeroVector() );

    // Initialize model output
    QUESO::GslVector model_output( full_model.get_observations() );
    model_output.cwSet(0.0);

    NitridationCalibration::FullModelLikelihood<QUESO::GslVector,QUESO::GslMatrix>
      likelihood( full_model );

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

    QUESO::GslVector step_size(full_model.get_model().param_space().zeroVector());
    unsigned int step_size_size = model_input.vector_variable_size("Optimizer/step_size");
    if( step_size_size != step_size.sizeLocal() )
      {
        if( env->fullRank() == 0 )
          {
            std::cerr << "Error: Initial step size mismatch!" << std::endl
                      << "input step size size = " << step_size_size << std::endl
                      << "sip step size size   = " << step_size.sizeLocal() << std::endl;
          }
        env.reset();
        MPI_Finalize();
        return 1;
       }

    for( unsigned int i = 0; i < guess_size; i++ )
       {
         step_size[i] = model_input("Optimizer/step_size", 0.0, i);
       }

     QUESO::OptimizerMonitor monitor(full_model.get_model().param_space().env(),10000);
     monitor.set_display_output(true,true);

     QUESO::GslOptimizer optimizer(likelihood);

     std::string solver_type = model_input("Optimizer/solver_type", "DIE");

     double h = model_input("Optimizer/finite_difference_step_size", 1.0e-8);
     optimizer.set_solver_type(solver_type);
     optimizer.set_step_size(step_size);
     optimizer.setFiniteDifferenceStepSize(h);


     if( env->fullRank() == 0 )
       {
         std::cout << std::endl
                   << "=============================================================" << std::endl
                   << "      Solving using: " << solver_type << std::endl
                   << "=============================================================" << std::endl;
       }


     optimizer.setInitialPoint(guess);

     optimizer.minimize(&monitor);

     const QUESO::GslVector& minimizer = optimizer.minimizer();

     for( unsigned int i = 0; i < guess_size; i++ )
       {
         if( env->fullRank() == 0 )
           {
             std::cout << "minimizer["<<i<<"] = " << minimizer[i] << std::endl;
           }
       }
  }

  MPI_Finalize();

#else

  std::cout << "Must have linked against a valid QUESO installation for this program to run." << std::endl;

#endif //NITCAL_HAVE_QUESO

  return 0;
}
