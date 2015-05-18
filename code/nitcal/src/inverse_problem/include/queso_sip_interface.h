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

    QuesoStatisticalInverseProblemInterface( const QUESO::BaseEnvironment& env,
                                             const std::string& method );

    virtual ~QuesoStatisticalInverseProblemInterface(){};

    void solve( const Vec& initial_guess );

    const QUESO::VectorSpace<Vec,Mat>&  get_param_space() const;
    const QUESO::VectorSubset<Vec,Mat>&  get_param_domain() const;
    const QUESO::BaseVectorRV<Vec,Mat>& get_prior_rv() const;
    const QUESO::GenericVectorRV<Vec,Mat>&  get_posterior_rv() const;
    const QUESO::BaseScalarFunction<Vec,Mat>& get_likelihood_func() const;

  protected:

    void create_proposal_cov_mat();

    void create_sip();

    void create_prior();
    void create_posterior();

    const QUESO::BaseEnvironment& _queso_env;

    std::string _method;

    QUESO::VectorSpace<Vec,Mat> *_param_space;
    QUESO::VectorSubset<Vec,Mat> *_param_domain;

    boost::scoped_ptr<QUESO::BaseScalarFunction<Vec,Mat> > _likelihood;

    boost::scoped_ptr<QUESO::BaseVectorRV<Vec,Mat> > _prior;
    boost::scoped_ptr<QUESO::GenericVectorRV<Vec,Mat> > _posterior;

    boost::scoped_ptr<QUESO::StatisticalInverseProblem<Vec,Mat> > _ip;

    boost::scoped_ptr<Mat> _proposal_cov_mat;

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

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO

#endif // NITCAL_QUESO_SIP_INTERFACE_H
