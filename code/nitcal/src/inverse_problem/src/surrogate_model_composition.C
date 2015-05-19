//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// This class
#include "surrogate_model_composition.h"

//NitCal
#include "model_builder.h"

// QUESO
#include <queso/InterpolationSurrogateIOASCII.h>
#include <queso/LinearLagrangeInterpolationSurrogate.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace NitridationCalibration
{
  template<class Vec,class Mat>
  SurrogateModelComposition<Vec,Mat>::SurrogateModelComposition(const QUESO::BaseEnvironment& queso_env,
                                                                const GetPot& model_input)
    : _model(ModelBuilder<Vec,Mat>::build_model(queso_env,model_input))
  {
    const QUESO::FullEnvironment& full_env = dynamic_cast<const QUESO::FullEnvironment&>( queso_env );

    // Grab the datasets we'll be working with
    unsigned int n_datasets = model_input.vector_variable_size("Likelihood/datasets");

    // Parse dataset names
    std::vector<std::string> dataset_names(n_datasets);
    for( unsigned int d = 0; d < n_datasets; d++ )
        dataset_names[d] = model_input( "Likelihood/datasets", "DIE!", d );

    // Setup forward model inputs
    _forward_model_inputs.resize(n_datasets,NULL);
    for( unsigned int d = 0; d < n_datasets; d++ )
      _forward_model_inputs[d] = new GetPot( dataset_names[d]+".in" );

    // Setup data space, 2 datapoints per dataset
    unsigned int n_datapoints = 2*n_datasets;
    QUESO::VectorSpace<Vec,Mat> data_space( queso_env, "data_", n_datapoints, NULL);

    _observations.reset( data_space.newVector() );
    _covariance.reset( data_space.newVector() );

    // Parse input files for observations, covariance
    for( unsigned int d = 0; d < n_datasets; d++ )
      {
        (*_observations)[2*d]   = (*_forward_model_inputs[d])("MassLossLikelihood/data_value", 0.0);
        (*_observations)[2*d+1] = (*_forward_model_inputs[d])("AverageNLikelihood/data_value", 0.0);

        double mass_loss_sigma = (*_forward_model_inputs[d])("MassLossLikelihood/sigma", 0.0);
        double avg_N_sigma = (*_forward_model_inputs[d])("AverageNLikelihood/sigma", 0.0);

        (*_covariance)[2*d]   = mass_loss_sigma*mass_loss_sigma;
        (*_covariance)[2*d+1] = avg_N_sigma*avg_N_sigma;
      }

    // Now setup IO classes. These will hold the InterpolationSurrodateData we need for evaluation
    _interp_io.resize( n_datapoints );
    for( unsigned int d = 0; d < n_datasets; d++ )
      {
        std::string mass_loss_filename = dataset_names[d]+"_massloss.dat";
        std::string avg_N_filename = dataset_names[d]+"_avgN.dat";

        _interp_io[2*d] = new QUESO::InterpolationSurrogateIOASCII<Vec,Mat>;
        _interp_io[2*d+1] = new QUESO::InterpolationSurrogateIOASCII<Vec,Mat>;

        _interp_io[2*d]->read( mass_loss_filename, full_env, "");
        _interp_io[2*d+1]->read( avg_N_filename, full_env, "");
      }

    // Now setup interpolation surrogates
    _interp_surrogate.resize( n_datapoints );
    for( unsigned int d = 0; d < n_datasets; d++ )
      {
        _interp_surrogate[2*d] =
          new QUESO::LinearLagrangeInterpolationSurrogate<Vec,Mat>( _interp_io[2*d]->data() );

        _interp_surrogate[2*d+1] =
          new QUESO::LinearLagrangeInterpolationSurrogate<Vec,Mat>( _interp_io[2*d+1]->data() );
      }
  }

  template<class Vec,class Mat>
  SurrogateModelComposition<Vec,Mat>::~SurrogateModelComposition()
  {
    for( std::vector<GetPot*>::iterator it = _forward_model_inputs.begin();
         it != _forward_model_inputs.end(); ++it )
      delete *it;

    for( typename std::vector<QUESO::InterpolationSurrogateIOASCII<Vec,Mat>*>::iterator it = _interp_io.begin();
         it != _interp_io.end(); ++it )
      delete *it;

    for( typename std::vector<QUESO::LinearLagrangeInterpolationSurrogate<Vec,Mat>*>::iterator it = _interp_surrogate.begin();
         it != _interp_surrogate.end(); ++it )
      delete *it;
  }

  template<class Vec,class Mat>
  void SurrogateModelComposition<Vec,Mat>::compute_values( const Vec& param_values,
                                                           Vec& model_output ) const
  {
    queso_require_equal_to( _observations->sizeGlobal(), model_output.sizeGlobal() );

    try
      {
        for( unsigned int d = 0; d < model_output.sizeGlobal(); d++ )
          model_output[d] = _interp_surrogate[d]->evaluate( param_values );
      }
    catch(...)
      {
        std::cout << "Caught exception in SurrogateModelComposition::compute_values()" << std::endl
                  << "param_values = " << param_values << std::endl
                  << "dim = " << _interp_io[0]->data().dim() << std::endl
                  << " x_min[0] = " << _interp_io[0]->data().x_min(0) << std::endl
                  << " x_max[0] = " << _interp_io[0]->data().x_max(0) << std::endl
                  << " x_min[1] = " << _interp_io[0]->data().x_min(1) << std::endl
                  << " x_max[1] = " << _interp_io[0]->data().x_max(1) << std::endl
                  << " spacing[0] = " << _interp_io[0]->data().spacing(0) << std::endl
                  << " spacing[1] = " << _interp_io[0]->data().spacing(1) << std::endl;

        MPI_Abort(MPI_COMM_WORLD,1);
      }
  }

  // Instantiate GSL version of this class
  template class SurrogateModelComposition<QUESO::GslVector,QUESO::GslMatrix>;

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
