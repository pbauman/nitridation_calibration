//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

#ifndef NITCAL_FULL_MODEL_LIKELIHOOD_H
#define NITCAL_FULL_MODEL_LIKELIHOOD_H

#include <queso_gaussian_likelihood_diagonal_covariance_interface.h>


namespace NitridationCalibration
{
  template<class Vec,class Mat>
  class FullModelComposition;

  template<class Vec,class Mat>
  class FullModelLikelihood :
    public QuesoGaussianLikelihoodDiagonalCovarianceInterface<Vec,Mat>
  {
  public:

    FullModelLikelihood( FullModelComposition<Vec,Mat>& full_model );

    virtual ~FullModelLikelihood(){};

    //! Override to compute the constants
    virtual void evaluateModel( const Vec& domainVector,
                                const Vec* domainDirection,
                                Vec& modelOutput,
                                Vec* gradVector,
                                Mat* hessianMatrix,
                                Vec* hessianEffect ) const;

  protected:

    FullModelComposition<Vec,Mat>& _full_model;

  };

} // end namespace NitridationCalibration

#endif // NITCAL_FULL_MODEL_LIKELIHOOD_H

#endif // NITCAL_HAVE_QUESO
