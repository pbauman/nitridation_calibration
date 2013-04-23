//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$ 
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef NITCAL_LIKELIHOOD_BASE_H
#define NITCAL_LIKELIHOOD_BASE_H

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// QUESO
#include "uqScalarFunction.h"

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  class LikelihoodBase : public uqBaseScalarFunctionClass<Vec,Mat>
  {
  public:

    LikelihoodBase( const char* prefix, 
		    const uqVectorSetClass<Vec,Mat>& domain_set,
		    const bool returns_ln);

    virtual ~LikelihoodBase();

    double actualValue(const Vec& domainVector, 
		       const Vec* domainDirection, 
		       Vec* gradVector, 
		       Mat* hessianMatrix, 
		       Vec* hessianEffect) const;

    double lnValue(const Vec& domainVector, 
		   const Vec* domainDirection, 
		   Vec* gradVector, 
		   Mat* hessianMatrix, 
		   Vec* hessianEffect) const;


    // This needs to be a const function since it sits in 
    // const functions.
    virtual double evaluate_likelihood() const = 0;

    virtual void update_parameters( const std::vector<double>& params ) const = 0;

  protected:

    // We inherit: m_env, m_prefix, m_domainSet

    // Users tells QUESO whether or not they're computing
    // log( likelihood_value )
    bool m_returns_ln;

    // Silly variable to monitor sampling
    // since this problem takes a bit.
    mutable unsigned int _sample_count;

  private:

    LikelihoodBase();

    void initialize();

  };

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO

#endif // NITCAL_LIKELIHOOD_BASE_H
