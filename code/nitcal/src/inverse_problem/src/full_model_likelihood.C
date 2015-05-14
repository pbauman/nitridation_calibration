//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// This class
#include "full_model_likelihood.h"

// NitCal
#include "full_model_composition.h"

// QUESO
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace NitridationCalibration
{
  template<class Vec,class Mat>
  FullModelLikelihood<Vec,Mat>::FullModelLikelihood( FullModelComposition<Vec,Mat>& full_model )
    : QuesoGaussianLikelihoodDiagonalCovarianceInterface<Vec,Mat>("likelihood_",
                                                                  full_model.get_model().param_domain(),
                                                                  full_model.get_observations(),
                                                                  full_model.get_covariance()),
    _full_model(full_model)
  {}

  template<class Vec,class Mat>
  void FullModelLikelihood<Vec,Mat>::evaluateModel( const Vec& domainVector,
                                                    const Vec* /*domainDirection*/,
                                                    Vec& modelOutput,
                                                    Vec* /*gradVector*/,
                                                    Mat* /*hessianMatrix*/,
                                                    Vec* /*hessianEffect*/ ) const
  {
    unsigned int n_params = domainVector.sizeGlobal();
    std::vector<double> param_values(n_params);

    for( unsigned int i = 0; i < n_params; i++ )
      param_values[i] = domainVector[i];

    unsigned int n_observations = modelOutput.sizeGlobal();
    std::vector<double> model_output(n_observations);

    _full_model.compute_values( param_values, model_output );

    for( unsigned int i = 0; i < n_observations; i++ )
      modelOutput[i] = model_output[i];
  }

  // Instantiate GSL version of this class
  template class FullModelLikelihood<QUESO::GslVector,QUESO::GslMatrix>;

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
