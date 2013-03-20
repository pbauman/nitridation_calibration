//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifdef HAVE_QUESO

// This class
#include "likelihood_base.h"

// QUESO
#include "uqGslVector.h"
#include "uqGslMatrix.h"

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  LikelihoodBase<Vec,Mat>::LikelihoodBase( const char* prefix, 
					   const uqVectorSetClass<Vec,Mat>& domain_set,
					   const bool returns_ln )
    : uqBaseScalarFunctionClass<Vec,Mat>( prefix, domain_set ),
      m_returns_ln( returns_ln ),
      _sample_count(0)
  {
    return;
  }

  template<class Vec,class Mat>
  LikelihoodBase<Vec,Mat>::~LikelihoodBase()
  {
    return;
  }


  template<class Vec,class Mat>
  double LikelihoodBase<Vec,Mat>::actualValue( const Vec& domainVector, 
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

    if(m_returns_ln) 
      {
#ifdef QUESO_EXPECTS_LN_LIKELIHOOD_INSTEAD_OF_MINUS_2_LN
	value = std::exp(value);
#else
	value = std::exp(-.5*value);
#endif
      }

    return value;
  }

  template<class Vec,class Mat>
  double LikelihoodBase<Vec,Mat>::lnValue( const Vec& domainVector, 
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

    // Pulled this code from uqGenericScalarFunctionClass<V,M>::lnValue
    if(m_returns_ln == false) 
      {

#ifdef QUESO_EXPECTS_LN_LIKELIHOOD_INSTEAD_OF_MINUS_2_LN
	value = std::log(value);
#else
	value = -2.*std::log(value);
#endif
      }

    return value;
  }

  // Instantiate GSL version of this class
  template class LikelihoodBase<uqGslVectorClass,uqGslMatrixClass>;

} // end namespace NitridationCalibration

#endif // HAVE_QUESO
