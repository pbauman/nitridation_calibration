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

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  QuesoStatisticalInverseProblemInterface<Vec,Mat>::QuesoStatisticalInverseProblemInterface( const std::string& method,
                                                                                             const QUESO::BaseEnvironment& queso_env,
                                                                                             const QUESO::BaseVectorRV<Vec,Mat>& prior,
                                                                                             const QUESO::BaseScalarFunction<Vec,Mat>& likelihood )
    : _method(method),
      _queso_env(queso_env),
      _posterior(NULL),
      _ip(NULL),
      _proposal_cov_mat(NULL)
  {
    // First create posterior
    this->_posterior.reset( new QUESO::GenericVectorRV<Vec,Mat>("post_", // Extra prefix before the default "rv_" prefix
                                                                prior.imageSet().vectorSpace() ) );

    // Now instantiate inverse problem
    _ip.reset( new QUESO::StatisticalInverseProblem<Vec,Mat>("", // No extra prefix before the default "ip_" prefix
                                                             NULL,
                                                             prior,
                                                             likelihood,
                                                             *_posterior) );

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
    if( _method == "metropolis-hastings")
      {
        _proposal_cov_mat.reset( _ip->priorRv().imageSet().vectorSpace().newDiagMatrix(1.0) );

        queso_require_equal_to( initial_guess.sizeGlobal(), _proposal_cov_mat->numRowsGlobal() );
        queso_require_equal_to( initial_guess.sizeGlobal(), _proposal_cov_mat->numCols() );

        for( unsigned int i = 0; i < initial_guess.sizeGlobal(); i++ )
          (*_proposal_cov_mat)(i,i) = initial_guess[i];

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

  /* ------------------------- Instantiate -------------------------*/
  template class QuesoStatisticalInverseProblemInterface<QUESO::GslVector,QUESO::GslMatrix>;

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
