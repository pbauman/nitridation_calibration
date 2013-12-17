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

int main(int argc, char* argv[])
{

  // Check command line count.
  if( argc < 3 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify QUESO input file." << std::endl;
      exit(1); // TODO: something more sophisticated for parallel runs?
    }

  std::string QUESO_input = argv[1];

  MPI_Init(&argc,&argv);

  uqFullEnvironmentClass* env = new uqFullEnvironmentClass(MPI_COMM_WORLD,
							   QUESO_input.c_str(),
							   "",
							   NULL );

  {
    NitridationCalibration::ConstantGammaNConstantGammaCNSIP<uqGslVectorClass,uqGslMatrixClass> sip( env, "multilevel",
                                                                                       argc, argv,
                                                                                       QUESO_input );

    const uqBaseScalarFunctionClass<uqGslVectorClass,uqGslMatrixClass>& raw_likelihood = sip.get_likelihood_func();

    const NitridationCalibration::ConstantGammaNConstantGammaCNLikelihood<uqGslVectorClass,uqGslMatrixClass>& likelihood = 
      libmesh_cast_ref<const NitridationCalibration::ConstantGammaNConstantGammaCNLikelihood<uqGslVectorClass,uqGslMatrixClass>& >( raw_likelihood );

  }
  //************************************************
  // Finalize environments
  //************************************************
  delete env;

  MPI_Finalize();

} // main

#else

int main()
{
  // If not built with QUESO, report that we skipped this test
  return 77;
}

#endif // NITCAL_HAVE_QUESO
