//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// This class
#include "surrogate_model_likelihood.h"

// NitCal
#include "surrogate_model_composition.h"

// QUESO
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace NitridationCalibration
{
  template<class Vec,class Mat>
  SurrogateModelLikelihood<Vec,Mat>::SurrogateModelLikelihood( SurrogateModelComposition<Vec,Mat>& surrogate_model )
    : QuesoGaussianLikelihoodDiagonalCovarianceInterface<Vec,Mat>("likelihood_",
                                                                  surrogate_model.get_model().param_domain(),
                                                                  surrogate_model.get_observations(),
                                                                  surrogate_model.get_covariance()),
    _surrogate_model(surrogate_model)
  {}

  template<class Vec,class Mat>
  void SurrogateModelLikelihood<Vec,Mat>::evaluateModel( const Vec& domainVector,
                                                         const Vec* /*domainDirection*/,
                                                         Vec& modelOutput,
                                                         Vec* /*gradVector*/,
                                                         Mat* /*hessianMatrix*/,
                                                         Vec* /*hessianEffect*/ ) const
  {
    _surrogate_model.compute_values( domainVector, modelOutput );
  }

  // Instantiate GSL version of this class
  template class SurrogateModelLikelihood<QUESO::GslVector,QUESO::GslMatrix>;

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
