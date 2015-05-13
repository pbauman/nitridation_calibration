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
#include <constant_gamma_n_constant_gamma_cn_model.h>
#include <arrhenius_gamma_n_constant_gamma_cn_model.h>
#include <arrhenius_gamma_n_arrhenius_gamma_cn_model.h>
#include <model_interpolation_builder.h>


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

  std::string forward_run_input = argv[1];
  std::string model_inputfile = argv[2];
  std::string QUESO_input = argv[3];

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

    boost::scoped_ptr<NitridationCalibration::ModelInterfaceBase<QUESO::GslVector,QUESO::GslMatrix> >
      model;

    GetPot model_input( model_inputfile );

    std::string model_type = model_input( "ModelType/model", "DIE!" );

    if( model_type == std::string("constant_gamma_cn") )
      model.reset( new NitridationCalibration::ConstantGammaCNModel<QUESO::GslVector,QUESO::GslMatrix>(argc,argv,*env,forward_run_input,model_inputfile ) );

      model.reset( new NitridationCalibration::ConstantGammaNConstantGammaCNModel<QUESO::GslVector,QUESO::GslMatrix>(argc,argv,*env,forward_run_input,model_inputfile ) );
    if( model_type == std::string("constant_gamma_n_constant_gamma_cn") )

    else if( model_type == std::string("arrhenius_gamma_n_constant_gamma_cn") )
      model.reset( new NitridationCalibration::ArrheniusGammaNConstantGammaCNModel<QUESO::GslVector,QUESO::GslMatrix>(argc,argv,*env,forward_run_input,model_inputfile ) );

    else if( model_type == std::string("arrhenius_gamma_n_constant_gamma_cn") )
      model.reset( new NitridationCalibration::ArrheniusGammaNArrheniusGammaCNModel<QUESO::GslVector,QUESO::GslMatrix>(argc,argv,*env,forward_run_input,model_inputfile ) );

    else
      {
        std::cerr << "Error: Invalid Model type! Found " << model_type << std::endl;
        env.reset();
        MPI_Finalize();
        return 1;
      }

    std::vector<unsigned int> n_points(model->param_domain().vectorSpace().dimGlobal());

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
      data(model->param_domain(),n_points,n_datasets);

    NitridationCalibration::ModelInterpolationBuilder<QUESO::GslVector,QUESO::GslMatrix>
      builder( data, *(model.get()) );

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
