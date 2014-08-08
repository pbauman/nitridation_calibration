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
#include "queso/StatisticalInverseProblem.h"

// NitCal
#include "likelihood_comm_handler.h"

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  class QuesoStatisticalInverseProblemInterface
  {
  public:

    QuesoStatisticalInverseProblemInterface( QUESO::BaseEnvironment* env,
                                             const std::string method );

    virtual ~QuesoStatisticalInverseProblemInterface();

    void solve();
    
    const QUESO::VectorSpace<Vec,Mat>&  get_param_space() const;
    const QUESO::VectorSubset<Vec,Mat>&  get_param_domain() const;
    const QUESO::BaseVectorRV<Vec,Mat>& get_prior_rv() const;
    const QUESO::GenericVectorRV<Vec,Mat>&  get_posterior_rv() const;
    const QUESO::BaseScalarFunction<Vec,Mat>& get_likelihood_func() const;

    const LikelihoodCommHandler& comm_handler() const;

  protected:

    void create_proposal_cov_mat();
  
    void create_sip();

    void create_prior();
    void create_posterior();

    QUESO::BaseEnvironment *_queso_env;

    QUESO::VectorSpace<Vec,Mat> *_param_space;
    QUESO::VectorSubset<Vec,Mat> *_param_domain;

    QUESO::BaseScalarFunction<Vec,Mat> *_likelihood;

    QUESO::BaseVectorRV<Vec,Mat> *_prior;
    QUESO::GenericVectorRV<Vec,Mat> *_posterior;
  
    QUESO::StatisticalInverseProblem<Vec,Mat>* _ip;

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
  const QUESO::VectorSpace<Vec,Mat>& QuesoStatisticalInverseProblemInterface<Vec,Mat>::get_param_space() const
  {
    return *_param_space;
  }

  template<class Vec,class Mat>
  inline
  const QUESO::VectorSubset<Vec,Mat>& QuesoStatisticalInverseProblemInterface<Vec,Mat>::get_param_domain() const
  {
    return *_param_domain;
  }

  template<class Vec,class Mat>
  inline
  const QUESO::BaseVectorRV<Vec,Mat>& QuesoStatisticalInverseProblemInterface<Vec,Mat>::get_prior_rv() const
  {
    return *_prior;
  }

  template<class Vec,class Mat>
  inline
  const QUESO::GenericVectorRV<Vec,Mat>& QuesoStatisticalInverseProblemInterface<Vec,Mat>::get_posterior_rv() const
  {
    return *_posterior;
  }

  template<class Vec,class Mat>
  inline
  const QUESO::BaseScalarFunction<Vec,Mat>& QuesoStatisticalInverseProblemInterface<Vec,Mat>::get_likelihood_func() const
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
