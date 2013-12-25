//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef NITCAL_QUESO_SIP_INTERFACE_H
#define NITCAL_QUESO_SIP_INTERFACE_H

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// Boost
#include "boost/scoped_ptr.hpp"

// Queso
#include "uqStatisticalInverseProblem.h"

// NitCal
#include "likelihood_comm_handler.h"

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  class QuesoStatisticalInverseProblemInterface
  {
  public:

    QuesoStatisticalInverseProblemInterface( uqBaseEnvironmentClass* env,
                                             const std::string method );

    virtual ~QuesoStatisticalInverseProblemInterface();

    void solve();
    
    const uqVectorSpaceClass<Vec,Mat>&  get_param_space() const;
    const uqVectorSubsetClass<Vec,Mat>&  get_param_domain() const;
    const uqBaseVectorRVClass<Vec,Mat>& get_prior_rv() const;
    const uqGenericVectorRVClass<Vec,Mat>&  get_posterior_rv() const;
    const uqBaseScalarFunctionClass<Vec,Mat>& get_likelihood_func() const;

    const LikelihoodCommHandler& comm_handler() const;

  protected:

    void create_proposal_cov_mat();
  
    void create_sip();

    void create_prior();
    void create_posterior();

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

    boost::scoped_ptr<LikelihoodCommHandler> _comm_handler;

  private:

    QuesoStatisticalInverseProblemInterface();

  };

  /* ------------------------- Inline Functions -------------------------*/

  template<class Vec,class Mat>
  inline
  const uqVectorSpaceClass<Vec,Mat>& QuesoStatisticalInverseProblemInterface<Vec,Mat>::get_param_space() const
  {
    return *_param_space;
  }

  template<class Vec,class Mat>
  inline
  const uqVectorSubsetClass<Vec,Mat>& QuesoStatisticalInverseProblemInterface<Vec,Mat>::get_param_domain() const
  {
    return *_param_domain;
  }

  template<class Vec,class Mat>
  inline
  const uqBaseVectorRVClass<Vec,Mat>& QuesoStatisticalInverseProblemInterface<Vec,Mat>::get_prior_rv() const
  {
    return *_prior;
  }

  template<class Vec,class Mat>
  inline
  const uqGenericVectorRVClass<Vec,Mat>& QuesoStatisticalInverseProblemInterface<Vec,Mat>::get_posterior_rv() const
  {
    return *_posterior;
  }

  template<class Vec,class Mat>
  inline
  const uqBaseScalarFunctionClass<Vec,Mat>& QuesoStatisticalInverseProblemInterface<Vec,Mat>::get_likelihood_func() const
  {
    return *_likelihood;
  }

  template<class Vec,class Mat>
  inline
  const LikelihoodCommHandler& QuesoStatisticalInverseProblemInterface<Vec,Mat>::comm_handler() const
  {
    return *(_comm_handler.get());
  }

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO

#endif // NITCAL_QUESO_SIP_INTERFACE_H
