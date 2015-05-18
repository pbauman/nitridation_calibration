//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// C++
#include <cmath>

// This class
#include "queso_sip_interface.h"

// QUESO
#include "queso/GslVector.h"
#include "queso/GslMatrix.h"
#include "queso/UniformVectorRV.h"

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  QuesoStatisticalInverseProblemInterface<Vec,Mat>::QuesoStatisticalInverseProblemInterface( const QUESO::BaseEnvironment& env,
                                                                                             const std::string& method )
    : _queso_env(env),
      _method(method),
      _param_space(NULL),
      _param_domain(NULL),
      _likelihood(NULL),
      _prior(NULL),
      _posterior(NULL),
      _ip(NULL),
      _proposal_cov_mat(NULL)
  {}

  template<class Vec,class Mat>
  void QuesoStatisticalInverseProblemInterface<Vec,Mat>::create_sip( )
  {
    _ip.reset( new QUESO::StatisticalInverseProblem<Vec,Mat>("", // No extra prefix before the default "ip_" prefix
                                                             NULL,
                                                             *_prior,
                                                             *_likelihood,
                                                             *_posterior) );

    return;
  }

  template<class Vec,class Mat>
  void QuesoStatisticalInverseProblemInterface<Vec,Mat>::solve( const Vec& initial_guess )
  {
    if (_queso_env.fullRank() == 0)
      {
	std::cout << "Beginning run of QUESO Statistical Inverse Problem using the "
		  << _method << " method."
		  << std::endl;
      }

    //******************************************************
    // Solve the inverse problem
    //******************************************************
    if( _method == "metropolis_hastings")
      {
	_ip->solveWithBayesMetropolisHastings( NULL, initial_guess, _proposal_cov_mat.get());
      }
    else if( _method == "multilevel" )
      {
	_ip->solveWithBayesMLSampling();
      }
    else
      {
	//FIXME: Proper stopping needed here.
	std::cout << "Invalid method " << _method << " for statistical inverse problem." << std::endl;
	MPI_Abort(MPI_COMM_WORLD,1);
      }

    if (_queso_env.fullRank() == 0)
      {
	std::cout << "Finishing run of QUESO Statistical Inverse Problem"
		  << std::endl;
      }

    return;
  }

  template<class Vec,class Mat>
  void QuesoStatisticalInverseProblemInterface<Vec,Mat>::create_prior()
  {
    this->_prior.reset( new QUESO::UniformVectorRV<Vec,Mat>("prior_", // Extra prefix before the default "rv_"
                                                            *(this->_param_domain) ) );

  }

  template<class Vec,class Mat>
  void QuesoStatisticalInverseProblemInterface<Vec,Mat>::create_posterior()
  {
    this->_posterior.reset( new QUESO::GenericVectorRV<Vec,Mat>("post_", // Extra prefix before the default "rv_" prefix
                                                                *(this->_param_space) ) );
  }

  /* ------------------------- Instantiate -------------------------*/
  template class QuesoStatisticalInverseProblemInterface<QUESO::GslVector,QUESO::GslMatrix>;

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
