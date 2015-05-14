//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

#include <queso_gaussian_likelihood_diagonal_covariance_interface.h>

#ifndef NITCAL_FULL_MODEL_LIKELIHOOD_H
#define NITCAL_FULL_MODEL_LIKELIHOOD_H

namespace NitridationCalibration
{
  template<class Vec,class Mat>
  class FullModelLikelihood :
    public QuesoGaussianLikelihoodDiagonalCovarianceInterface<Vec,Mat>
  {
  public:

    FullModelLikelihood( FullModelComposition<Vec,Mat>& full_model );

    virtual ~FullModelLikelihood(){};

    //! Override to compute the constants
    virtual double evaluate_model( const Vec& domainVector,
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
