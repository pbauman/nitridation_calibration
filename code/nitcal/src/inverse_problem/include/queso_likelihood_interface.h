//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$ 
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef NITCAL_QUESO_LIKELIHOOD_INTERFACE_H
#define NITCAL_QUESO_LIKELIHOOD_INTERFACE_H

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// QUESO
#include "queso/ScalarFunction.h"

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  class QuesoLikelihoodInterface : public QUESO::BaseScalarFunction<Vec,Mat>
  {
  public:

    QuesoLikelihoodInterface( const char* prefix, 
                              const QUESO::VectorSet<Vec,Mat>& domain_set,
                              const bool returns_ln);

    virtual ~QuesoLikelihoodInterface();

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

  private:

    QuesoLikelihoodInterface();

    void initialize();

  };

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO

#endif // NITCAL_QUESO_LIKELIHOOD_INTERFACE_H
