//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// This class
#include "full_model_composition.h"

//NitCal
#include "full_model_evaluator.h"
#include "model_builder.h"

// QUESO
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace NitridationCalibration
{
  template<class Vec,class Mat>
  FullModelComposition<Vec,Mat>::FullModelComposition( int argc, char** argv,
                                                       const QUESO::BaseEnvironment& queso_env,
                                                       const GetPot& model_input )
    : _model(ModelBuilder<Vec,Mat>::build_model(queso_env,model_input)),
      _comm_handler(queso_env.subComm().Comm(),
                    model_input.vector_variable_size("Likelihood/datasets") )
  {
    // Grab the datasets we'll be working with
    unsigned int n_datasets = model_input.vector_variable_size("Likelihood/datasets");

    std::vector<std::string> datasets(n_datasets);
    for( unsigned int d = 0; d < n_datasets; d++ )
      {
        datasets[d] = model_input( "Likelihood/datasets", "DIE!", d );
      }

    // This is the dataset the current set of processors is going to work on
    int dataset_index = this->_comm_handler.get_dataset_index();

    // Input for this dataset
    _forward_run_input.reset( new GetPot(datasets[dataset_index]) );


    // Setup data space, 2 datapoints per dataset
    unsigned int n_datapoints = 2*n_datasets;
    QUESO::VectorSpace<Vec,Mat> data_space( queso_env, "data_", n_datapoints, NULL);

    _observations.reset( data_space.newVector() );
    _covariance.reset( data_space.newVector() );

    // Now parse data values and the corresponding covariances
    // Each processor parses its own dataset
    // Then we'll gather/broadcast to everyone
    std::vector<double> local_values(2);
    std::vector<double> all_values(n_datapoints);

    // Convention, mass_loss is first, then avg_N
    local_values[0] = (*_forward_run_input)("MassLossLikelihood/data_value", 0.0);
    local_values[1] = (*_forward_run_input)("AverageNLikelihood/data_value", 0.0);

    if( _comm_handler.get_inter0_rank() >= 0 )
      MPI_Gather( &local_values[0], 2, MPI_DOUBLE,
                  &all_values[0], 2, MPI_DOUBLE, 0,
                  _comm_handler.get_inter_chain_0_comm() );

    MPI_Bcast( &all_values[0], n_datapoints, MPI_DOUBLE,
               0, _comm_handler.get_inter_chain_comm() );

    for( unsigned int i = 0; i < n_datapoints; i++ )
      (*_observations)[i] = all_values[i];

    local_values[0] = (*_forward_run_input)("MassLossLikelihood/sigma", -1.0);
    local_values[1] = (*_forward_run_input)("AverageNLikelihood/sigma", -1.0);

    if( _comm_handler.get_inter0_rank() >= 0 )
      MPI_Gather( &local_values[0], 2, MPI_DOUBLE,
                  &all_values[0], 2, MPI_DOUBLE, 0,
                  _comm_handler.get_inter_chain_0_comm() );

    MPI_Bcast( &all_values[0], n_datapoints, MPI_DOUBLE,
               0, _comm_handler.get_inter_chain_comm() );

    for( unsigned int i = 0; i < n_datapoints; i++ )
      (*_covariance)[i] = all_values[i];


    // Now setup model to be evaluated on this set of processors
    // We do this last because of the UFO check in GRINS
    _model_evaluator.reset( new FullModelEvaluator<Vec,Mat>(argc,argv,
                                                            queso_env,
                                                            *(_forward_run_input.get()),
                                                            _comm_handler.get_split_chain_comm(),
                                                            *(_model.get())) );
  }

  template<class Vec,class Mat>
  void FullModelComposition<Vec,Mat>::compute_values( const std::vector<double>& param_values,
                                                      std::vector<double>& model_output ) const
  {
    // Each data set will compute the two outputs
    std::vector<double> local_values(2);
    _model_evaluator->compute_values( param_values, local_values );

    // Now gather the local_values from all evaluated datasets
    // to be able to insert into the model_output vector.
    // First, we gather from processor 0 from each set of workers
    queso_assert_equal_to( model_output.size(), _observations->sizeGlobal() );

    // We can only call this if we're a member of inter_chain_0
    // By convention, inter0_rank is negative if this processor
    // is not in that communicator
    if( _comm_handler.get_inter0_rank() >= 0 )
      MPI_Gather( &local_values[0], 2, MPI_DOUBLE,
                  &model_output[0], 2, MPI_DOUBLE, 0,
                  _comm_handler.get_inter_chain_0_comm() );

    // Now broadcast to the rest of the workers on the whole chain
    // We did the split, so processor 0 on inter_chain_0 is also
    // processor 0 on inter_chain
    MPI_Bcast( &model_output[0], _observations->sizeGlobal(), MPI_DOUBLE,
               0, _comm_handler.get_inter_chain_comm() );
  }

  // Instantiate GSL version of this class
  template class FullModelComposition<QUESO::GslVector,QUESO::GslMatrix>;

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
