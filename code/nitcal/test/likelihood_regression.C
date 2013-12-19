//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// C++
#include <string>
#include <iostream>

// NitCal
#include "constant_gamma_n_constant_gamma_cn_sip.h"
#include "constant_gamma_n_constant_gamma_cn_likelihood.h"

// QUESO
#include "uqEnvironment.h"
#include "uqGslVector.h"
#include "uqGslMatrix.h"

// GRINS
#include "grins/math_constants.h"

int main(int argc, char* argv[])
{

  // Check command line count.
  if( argc < 4 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify QUESO input file, mass loss value, and average N value" << std::endl;
      exit(1); // TODO: something more sophisticated for parallel runs?
    }

  std::string QUESO_input = argv[1];

  MPI_Init(&argc,&argv);

  uqFullEnvironmentClass* env = new uqFullEnvironmentClass(MPI_COMM_WORLD,
							   QUESO_input.c_str(),
							   "",
							   NULL );

  int return_flag = 0;

  {
    NitridationCalibration::ConstantGammaNConstantGammaCNSIP<uqGslVectorClass,uqGslMatrixClass> sip( env, "multilevel",
                                                                                       argc, argv,
                                                                                       QUESO_input );

    const uqBaseScalarFunctionClass<uqGslVectorClass,uqGslMatrixClass>& raw_likelihood = sip.get_likelihood_func();

    const NitridationCalibration::ConstantGammaNConstantGammaCNLikelihood<uqGslVectorClass,uqGslMatrixClass>& likelihood = 
      libmesh_cast_ref<const NitridationCalibration::ConstantGammaNConstantGammaCNLikelihood<uqGslVectorClass,uqGslMatrixClass>& >( raw_likelihood );


    const double computed_likelihood = likelihood.evaluate_likelihood();

    //std::cout << std::scientific << std::setprecision(16) << "likelihood = " << computed_likelihood << std::endl;

    const GetPot& sip_input = sip.get_sip_input();

    const GetPot& fp_input = sip.get_forward_run_input();

    const double mass_loss_data = fp_input( "MassLossLikelihood/data_value", 0.0 );
    if( !fp_input.have_variable( "MassLossLikelihood/data_value" ) )
      {
        std::cerr << "Error: Missing mass loss data_value in fp_input" << std::endl;
        return 1;
      }

    const double mass_loss_sigma = fp_input( "MassLossLikelihood/sigma", 0.0 );
    if( !fp_input.have_variable( "MassLossLikelihood/sigma" ) )
      {
        std::cerr << "Error: Missing mass loss sigma in fp_input" << std::endl;
        return 1;
      }

    const double avg_N_data = fp_input( "AverageNLikelihood/data_value", 0.0 );
    if( !fp_input.have_variable( "AverageNLikelihood/data_value" ) )
      {
        std::cerr << "Error: Missing average N data_value in fp_input" << std::endl;
        return 1;
      }

    const double avg_N_sigma = fp_input( "AverageNLikelihood/sigma", 0.0 );
    if( !fp_input.have_variable( "AverageNLikelihood/sigma" ) )
      {
        std::cerr << "Error: Missing average N sigma in fp_input" << std::endl;
        return 1;
      }

    const double mass_loss_value = atof(argv[2]);
    const double avg_N_value = atof(argv[3]);

    const double exact_ln_likelihood = std::log( 1.0/(mass_loss_sigma*std::sqrt(2*GRINS::Constants::pi) ) )
      + -0.5*(mass_loss_value - mass_loss_data)*(mass_loss_value - mass_loss_data)/(mass_loss_sigma*mass_loss_sigma)
      + std::log( 1.0/(avg_N_sigma*std::sqrt(2*GRINS::Constants::pi)) )
      + -0.5*(avg_N_value - avg_N_data)*(avg_N_value - avg_N_data)/(avg_N_sigma*avg_N_sigma);

    //std::cout << "exact_ln_likelihood = " << exact_ln_likelihood << std::endl;

    const double tol = 1.0e-10;

    if( std::fabs( (computed_likelihood - exact_ln_likelihood)/exact_ln_likelihood ) > tol )
      {
        std::cerr << std::scientific << std::setprecision(16)
                  << "Error: Tolerance exceeded in computed likelihood" << std::endl
                  << "       tolerance              = " << tol << std::endl
                  << "       computed_ln_likelihood = " << computed_likelihood << std::endl
                  << "       exact_ln_likelihood    = " << exact_ln_likelihood << std::endl;

        return_flag = 1;
      }

  }
  //************************************************
  // Finalize environments
  //************************************************
  delete env;

  MPI_Finalize();

  return return_flag;

} // main

#else

int main()
{
  // If not built with QUESO, report that we skipped this test
  return 77;
}

#endif // NITCAL_HAVE_QUESO
