//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// This class
#include "queso_likelihood_interface.h"

// QUESO
#include "queso/GslVector.h"
#include "queso/GslMatrix.h"

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  QuesoLikelihoodInterface<Vec,Mat>::QuesoLikelihoodInterface( const char* prefix, 
                                                               const QUESO::VectorSet<Vec,Mat>& domain_set,
                                                               const bool returns_ln )
    : QUESO::BaseScalarFunction<Vec,Mat>( prefix, domain_set ),
      m_returns_ln( returns_ln )
  {
    return;
  }

  template<class Vec,class Mat>
  QuesoLikelihoodInterface<Vec,Mat>::~QuesoLikelihoodInterface()
  {
    return;
  }


  template<class Vec,class Mat>
  double QuesoLikelihoodInterface<Vec,Mat>::actualValue( const Vec& domainVector, 
                                                         const Vec* /*domainDirection*/, 
                                                         Vec* /*gradVector*/, 
                                                         Mat* /*hessianMatrix*/, 
                                                         Vec* /*hessianEffect*/ ) const
  {
    //********************************************************************
    // Copy contents of domainVector to std::vector.
    // Doing this copy so we don't have to template the parameters class.
    //********************************************************************
    unsigned int num_params = domainVector.sizeLocal();

    std::vector<double> param_values( num_params );
  
    for(unsigned int i = 0; i < num_params; i++ )
      {
	param_values[i] = domainVector[i];
      }

    this->update_parameters( param_values );

    //********************************************************************
    // Evaluate likelihood
    //********************************************************************
    double value = this->evaluate_likelihood( );

    // If evaluate_likelihood is computing the ln, then we need to exp
    if(m_returns_ln)
      {
	value = std::exp(value);
      }

    return value;
  }

  template<class Vec,class Mat>
  double QuesoLikelihoodInterface<Vec,Mat>::lnValue( const Vec& domainVector, 
                                                     const Vec* /*domainDirection*/, 
                                                     Vec* /*gradVector*/, 
                                                     Mat* /*hessianMatrix*/, 
                                                     Vec* /*hessianEffect*/  ) const
  {
    //********************************************************************
    // Copy contents of domainVector to std::vector.
    // Doing this copy so we don't have to template the parameters class.
    //********************************************************************
    unsigned int num_params = domainVector.sizeLocal();

    std::vector<double> param_values( num_params );
  
    for(unsigned int i = 0; i < num_params; i++ )
      {
	param_values[i] = domainVector[i];
      }

    this->update_parameters( param_values );

    //********************************************************************
    // Evaluate likelihood
    //********************************************************************
    double value = this->evaluate_likelihood( );

    // // If evaluate_likelihood is *not* computing the ln, then we need to log
    if(m_returns_ln == false)
      {
	value = std::log(value);
      }

    return value;
  }

  // Instantiate GSL version of this class
  template class QuesoLikelihoodInterface<QUESO::GslVector,QUESO::GslMatrix>;

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
