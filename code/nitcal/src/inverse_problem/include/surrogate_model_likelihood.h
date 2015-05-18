//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

#ifndef NITCAL_SURROGATE_MODEL_LIKELIHOOD_H
#define NITCAL_SURROGATE_MODEL_LIKELIHOOD_H

#include <queso_gaussian_likelihood_diagonal_covariance_interface.h>


namespace NitridationCalibration
{
  template<class Vec,class Mat>
  class SurrogateModelComposition;

  template<class Vec,class Mat>
  class SurrogateModelLikelihood :
    public QuesoGaussianLikelihoodDiagonalCovarianceInterface<Vec,Mat>
  {
  public:

    SurrogateModelLikelihood( SurrogateModelComposition<Vec,Mat>& full_model );

    virtual ~SurrogateModelLikelihood(){};

    //! Override to compute the constants
    virtual void evaluateModel( const Vec& domainVector,
                                const Vec* domainDirection,
                                Vec& modelOutput,
                                Vec* gradVector,
                                Mat* hessianMatrix,
                                Vec* hessianEffect ) const;

  protected:

    SurrogateModelComposition<Vec,Mat>& _surrogate_model;

  };

} // end namespace NitridationCalibration

#endif // NITCAL_SURROGATE_MODEL_LIKELIHOOD_H

#endif // NITCAL_HAVE_QUESO
