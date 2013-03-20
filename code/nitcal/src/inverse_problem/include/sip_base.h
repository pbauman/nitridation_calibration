//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef NITCAL_SIP_BASE_H
#define NITCAL_SIP_BASE_H

#ifdef HAVE_QUESO

#include "uqStatisticalInverseProblem.h"

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  class StatisticalInverseProblemBase
  {
  public:

    StatisticalInverseProblemBase( uqBaseEnvironmentClass* env,
				   const std::string method );

    virtual ~StatisticalInverseProblemBase();

    void solve();
    
    const uqVectorSpaceClass<Vec,Mat>&  get_param_space() const;
    const uqVectorSubsetClass<Vec,Mat>&  get_param_domain() const;
    const uqBaseVectorRVClass<Vec,Mat>& get_prior_rv() const;
    const uqGenericVectorRVClass<Vec,Mat>&  get_posterior_rv() const;
    const uqBaseScalarFunctionClass<Vec,Mat>& get_likelihood_func() const;

  protected:

    void create_proposal_cov_mat();
  
    void create_sip();

    uqBaseEnvironmentClass *_queso_env;

    uqVectorSpaceClass<Vec,Mat> *_param_space;
    uqVectorSubsetClass<Vec,Mat> *_param_domain;

    uqBaseScalarFunctionClass<Vec,Mat> *_likelihood;

    uqBaseVectorRVClass<Vec,Mat> *_prior;
    uqGenericVectorRVClass<Vec,Mat> *_posterior;
  
    uqStatisticalInverseProblemClass<Vec,Mat>* _ip;

    Mat* _proposal_cov_mat;
    Vec* _param_initials;

    std::string _method;

  private:

    StatisticalInverseProblemBase();

  };

  /* ------------------------- Inline Functions -------------------------*/

  template<class Vec,class Mat>
  inline
  const uqVectorSpaceClass<Vec,Mat>& StatisticalInverseProblemBase<Vec,Mat>::get_param_space() const
  {
    return *_param_space;
  }

  template<class Vec,class Mat>
  inline
  const uqVectorSubsetClass<Vec,Mat>& StatisticalInverseProblemBase<Vec,Mat>::get_param_domain() const
  {
    return *_param_domain;
  }

  template<class Vec,class Mat>
  inline
  const uqBaseVectorRVClass<Vec,Mat>& StatisticalInverseProblemBase<Vec,Mat>::get_prior_rv() const
  {
    return *_prior;
  }

  template<class Vec,class Mat>
  inline
  const uqGenericVectorRVClass<Vec,Mat>& StatisticalInverseProblemBase<Vec,Mat>::get_posterior_rv() const
  {
    return *_posterior;
  }

  template<class Vec,class Mat>
  inline
  const uqBaseScalarFunctionClass<Vec,Mat>& StatisticalInverseProblemBase<Vec,Mat>::get_likelihood_func() const
  {
    return *_likelihood;
  }

} // end namespace NitridationCalibration

#endif // HAVE_QUESO

#endif // NITCAL_SIP_BASE_H
