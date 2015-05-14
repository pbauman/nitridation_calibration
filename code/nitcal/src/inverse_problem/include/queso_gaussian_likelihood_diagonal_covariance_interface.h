//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

#ifndef NITCAL_QUESO_GAUSSIAN_LIKELIHOOD_DIAGONAL_COVARIANCE_INTERFACE_H
#define NITCAL_QUESO_GAUSSIAN_LIKELIHOOD_DIAGONAL_COVARIANCE_INTERFACE_H

#include <queso/GaussianLikelihoodDiagonalCovariance.h>

namespace NitridationCalibration
{
  template<class Vec,class Mat>
   class QuesoGaussianLikelihoodDiagonalCovarianceInterface :
    public QUESO::GaussianLikelihoodDiagonalCovariance<Vec,Mat>
  {
  public:

    QuesoGaussianLikelihoodDiagonalCovarianceInterface( const char *prefix,
                                                        const QUESO::VectorSet<Vec,Mat>& domainSet,
                                                        const Vec& observations,
                                                        const Vec& covariance );

    virtual ~QuesoGaussianLikelihoodDiagonalCovarianceInterface(){};

    //! Override to compute the constants
    virtual double lnValue( const Vec& domainVector,
                            const Vec* domainDirection,
                            Vec* gradVector,
                            Mat* hessianMatrix,
                            Vec* hessianEffect ) const;

  private:

    QuesoGaussianLikelihoodDiagonalCovarianceInterface();

  };

} // end namespace NitridationCalibration

#endif // NITCAL_QUESO_GAUSSIAN_LIKELIHOOD_DIAGONAL_COVARIANCE_INTERFACE_H

#endif // NITCAL_HAVE_QUESO
