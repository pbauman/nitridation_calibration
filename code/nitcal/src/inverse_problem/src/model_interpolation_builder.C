//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// This class
#include <model_interpolation_builder.h>

// QUESO
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace NitridationCalibration
{
  template<class Vec,class Mat>
  ModelInterpolationBuilder<Vec,Mat>::ModelInterpolationBuilder( QUESO::InterpolationSurrogateDataSet<Vec,Mat>& data,
                                                                 const ModelEvaluatorBase<Vec,Mat>& model )
    : QUESO::InterpolationSurrogateBuilder<Vec,Mat>(data),
    _model(model),
    _count(0)
  {}

  template<class Vec,class Mat>
  void ModelInterpolationBuilder<Vec,Mat>::evaluate_model( const Vec& domainVector,
                                                           std::vector<double>& values )
  {
    // Only print out progress indicator on rank 0
    if( this->_model.get_model().param_space().env().fullRank() == 0 )
      std::cout << "ModelInterpolationBuilder evaluation " << _count << std::endl;

    //********************************************************************
    // Copy contents of domainVector to std::vector.
    // Doing this copy so we don't have to template the parameters class.
    //********************************************************************
    unsigned int num_params = domainVector.sizeGlobal();

    std::vector<double> param_values( num_params );
    for(unsigned int i = 0; i < num_params; i++ )
      {
	param_values[i] = domainVector[i];
      }

    _model.compute_values( param_values, values );

    _count++;
  }

  // Instantiate GSL version of this class
  template class ModelInterpolationBuilder<QUESO::GslVector,QUESO::GslMatrix>;

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
