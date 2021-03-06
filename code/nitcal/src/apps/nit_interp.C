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
#include <boost/scoped_ptr.hpp>

// NitCal
#include <model_interpolation_builder.h>
#include <full_model_composition.h>


#ifdef NITCAL_HAVE_QUESO
// QUESO
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/BoxSubset.h>
#include <queso/InterpolationSurrogateBuilder.h>
#include <queso/InterpolationSurrogateIOASCII.h>
#endif // HAVE_QUESO


int main(int argc, char* argv[])
{
#ifdef NITCAL_HAVE_QUESO

  // Check command line count.
  if( argc < 4 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify foward_run, model, and QUESO input files."
                << std::endl;
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

    std::vector<unsigned int> n_points(full_model.get_model().param_domain().vectorSpace().dimGlobal());

    unsigned int num_points = model_input.vector_variable_size( "ModelInterpolation/n_points");

    if( n_points.size() != num_points )
      {
        std::cerr << "Error: Mismatch in n_points and parameter space dim!"
                  << std::endl
                  << "num_points = " << num_points << std::endl
                  << "dim        = " << n_points.size()
                  << std::endl;
        env.reset();
        MPI_Finalize();
        return 1;
      }

    for( unsigned int n = 0; n < num_points; n++ )
      n_points[n] = model_input( "ModelInterpolation/n_points", 0, n );

    // Always two datasets: mass_loss, avg_N
    unsigned int n_datasets = 2;
    QUESO::InterpolationSurrogateDataSet<QUESO::GslVector, QUESO::GslMatrix>
      data(full_model.get_model().param_domain(),
           n_points,
           n_datasets);

    NitridationCalibration::ModelInterpolationBuilder<QUESO::GslVector,QUESO::GslMatrix>
      builder( data, full_model.get_model_evaluator() );

    // The expensive part. The builder will now evaluate the model for all the
    // desired points in parameter space. This will build both interpolants.
    builder.build_values();

    // Now that we've built the data, we write it out so we can reuse it later
    QUESO::InterpolationSurrogateIOASCII<QUESO::GslVector, QUESO::GslMatrix>
      data_writer;

    data_writer.write( "constant_gamma_cn_mass_loss.dat", data.get_dataset(0) );
    data_writer.write( "constant_gamma_cn_avg_N.dat", data.get_dataset(1) );
  }

  MPI_Finalize();

#else

  std::cout << "Must have linked against a valid QUESO installation for this program to run." << std::endl;

#endif //NITCAL_HAVE_QUESO

  return 0;
}
