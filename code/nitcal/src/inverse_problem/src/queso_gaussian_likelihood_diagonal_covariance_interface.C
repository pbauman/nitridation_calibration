//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

#include "queso_gaussian_likelihood_diagonal_covariance_interface.h"

// QUESO
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  QuesoGaussianLikelihoodDiagonalCovarianceInterface<Vec,Mat>::QuesoGaussianLikelihoodDiagonalCovarianceInterface( const char *prefix,
                                                                                                                   const QUESO::VectorSet<Vec,Mat>& domainSet,
                                                                                                                   const Vec& observations,
                                                                                                                   const Vec& covariance)
    : QUESO::GaussianLikelihoodDiagonalCovariance<Vec,Mat>(prefix,domainSet,observations,covariance)
  {}


  template<class Vec,class Mat>
  double QuesoGaussianLikelihoodDiagonalCovarianceInterface<Vec,Mat>::
  lnValue( const Vec& domainVector,
           const Vec* domainDirection,
           Vec* gradVector,
           Mat* hessianMatrix,
           Vec* hessianEffect ) const
  {
    /* Multivariate Gaussian PDF:
       (2 \pi)^{-k/2} \det{\Sigma}^{-1/2}
       \exp{ -1/2 (x-\mu)^T \Sigma^{-1} x(-\mu) }

       ==> ln(pdf) = -k/2*ln(2 \pi) -1/2*\det{\Simga}
                           - 1/2((x-\mu)^T \Sigma^{-1} x(-\mu))

    The base class computes the third term, so we need to compute the
    constants. */

    // -1/2 || f(m) - d ||^2
    double ln_value = QUESO::GaussianLikelihoodDiagonalCovariance<Vec,Mat>::lnValue(domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect);

    // Get the dimension and convert to double
    //double dim = this->m_observations.sizeGlobal();

    //const double pi = 3.1415926535897932384626433832795029L;
    //const double log_2_pi = std::log(2.0*pi);

    // Leading pi constant
    //ln_value -= dim/2.0*log_2_pi;

    // Since we have a diagonal covariance, the \det is just the
    // product of the covariances
    /*
    double det_sigma = 1.0;
    for( unsigned int i = 0; i < dim; i++ )
      {
        queso_assert_greater( this->m_covariance[i], 0.0 );
        det_sigma *= this->m_covariance[i];
      }
    */
    //ln_value -= 0.5*std::log(det_sigma);

    return ln_value;
  }

  // Instantiate GSL version of this class
  template class QuesoGaussianLikelihoodDiagonalCovarianceInterface<QUESO::GslVector,QUESO::GslMatrix>;

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
