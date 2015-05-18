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
#include "queso/GenericVectorRV.h"

// NitCal
#include "likelihood_comm_handler.h"

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  class QuesoStatisticalInverseProblemInterface
  {
  public:

    QuesoStatisticalInverseProblemInterface( const std::string& method,
                                             const QUESO::BaseEnvironment& queso_env,
                                             const QUESO::BaseVectorRV<Vec,Mat>& prior,
                                             const QUESO::BaseScalarFunction<Vec,Mat>& likelihood );

    virtual ~QuesoStatisticalInverseProblemInterface(){};

    void solve( const Vec& initial_guess );

    const QUESO::GenericVectorRV<Vec,Mat>&  get_posterior_rv() const;

  protected:

    std::string _method;

    const QUESO::BaseEnvironment& _queso_env;

    boost::scoped_ptr<QUESO::GenericVectorRV<Vec,Mat> > _posterior;

    boost::scoped_ptr<QUESO::StatisticalInverseProblem<Vec,Mat> > _ip;

    boost::scoped_ptr<Mat> _proposal_cov_mat;

  private:

    QuesoStatisticalInverseProblemInterface();

  };

  /* ------------------------- Inline Functions -------------------------*/
  template<class Vec,class Mat>
  inline
  const QUESO::GenericVectorRV<Vec,Mat>& QuesoStatisticalInverseProblemInterface<Vec,Mat>::get_posterior_rv() const
  {
    return *_posterior;
  }

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO

#endif // NITCAL_QUESO_SIP_INTERFACE_H
