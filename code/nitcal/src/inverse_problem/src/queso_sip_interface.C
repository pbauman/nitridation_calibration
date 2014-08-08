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
  QuesoStatisticalInverseProblemInterface<Vec,Mat>::QuesoStatisticalInverseProblemInterface( QUESO::BaseEnvironment *env,
                                                                                             const std::string method )
    : _queso_env(env),
      _method(method)
  {
    return;
  }

  template<class Vec,class Mat>
  QuesoStatisticalInverseProblemInterface<Vec,Mat>::~QuesoStatisticalInverseProblemInterface( )
  {
    delete _ip;
    
    return;
  }

  template<class Vec,class Mat>
  void QuesoStatisticalInverseProblemInterface<Vec,Mat>::create_sip( )
  {
    _ip = new QUESO::StatisticalInverseProblem<Vec,Mat>("", // No extra prefix before the default "ip_" prefix
							NULL,
							*_prior,
							*_likelihood,
							*_posterior);
    
    return;
  }
  
  template<class Vec,class Mat>
  void QuesoStatisticalInverseProblemInterface<Vec,Mat>::solve( )
  {
    if (_queso_env->fullRank() == 0) 
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
	// FIXME: Should add options to control with of the methods
	// FIXME: we use that are available in QUESO.
	_ip->solveWithBayesMetropolisHastings( NULL, *_param_initials, _proposal_cov_mat);
      }
    else if( _method == "multilevel" )
      {
	_ip->solveWithBayesMLSampling();
      }
    else
      {
	//FIXME: Proper stopping needed here.
	std::cout << "Invalid method " << _method << " for statistical inverse problem." << std::endl;
	exit(1);
      }
    
    /*
      
      
    //******************************************************
    // Write data to disk, to be used by 'sip_plot.m' afterwards
    //******************************************************
    if (_queso_env->fullRank() == 0) {
    std::cout << "Inverse problem solved. Writing data to disk now ...\n"
    << std::endl;
    }
    
    char varPrefixName[64+1];
    std::set<unsigned int> auxSet;
    auxSet.insert(0);
    
    sprintf(varPrefixName,"sip_appl_paramMeans");
    paramMeans.subWriteContents(varPrefixName,
    "outputData/appl_output",
    auxSet);
    sprintf(varPrefixName,"sip_appl_covMatrix");
    covMatrix->subWriteContents(varPrefixName,
    "outputData/appl_output",
    auxSet);
    sprintf(varPrefixName,"sip_appl_covMatrixInverse");
    covMatrixInverse->subWriteContents(varPrefixName,
    "outputData/appl_output",
    auxSet);

    //******************************************************
    // Write weighted squared norm to disk, to be used by 'sip_plot.m' afterwards
    //******************************************************
    // Define auxVec
    const QUESO::BaseVectorRealizer<Vec,Mat>& postRealizer = postRv.realizer();
    QUESO::VectorSpace<Vec,Mat> auxSpace(env,"",postRealizer.subPeriod(),NULL);
    Vec auxVec(auxSpace.zeroVector());

    // Populate auxVec
    Vec tmpVec (paramSpace.zeroVector());
    Vec diffVec(paramSpace.zeroVector());
    for (unsigned int i = 0; i < auxSpace.dimLocal(); ++i) {
    postRealizer.realization(tmpVec);
    diffVec = tmpVec - paramMeans;
    auxVec[i] = scalarProduct(diffVec, *covMatrixInverse * diffVec);
    }

    // Write auxVec to disk
    sprintf(varPrefixName,"sip_appl_d");
    auxVec.subWriteContents(varPrefixName,
    "outputData/appl_output",
    auxSet);

    //******************************************************
    // Release memory before leaving routine.
    //******************************************************
    delete covMatrixInverse;
    delete covMatrix;

    */

    if (_queso_env->fullRank() == 0) 
      {
	std::cout << "Finishing run of QUESO Statistical Inverse Problem"
		  << std::endl;
      }

    return;
  }

  template<class Vec,class Mat>
  void QuesoStatisticalInverseProblemInterface<Vec,Mat>::create_prior()
  {
    this->_prior =
      new QUESO::UniformVectorRV<Vec,Mat>("prior_", // Extra prefix before the default "rv_"
					  *(this->_param_domain) );

    return;
  }

  template<class Vec,class Mat>
  void QuesoStatisticalInverseProblemInterface<Vec,Mat>::create_posterior()
  {
    this->_posterior =
      new QUESO::GenericVectorRV<Vec,Mat>("post_", // Extra prefix before the default "rv_" prefix
					  *(this->_param_space) );

    return;
  }

  /* ------------------------- Instantiate -------------------------*/
  template class QuesoStatisticalInverseProblemInterface<QUESO::GslVector,QUESO::GslMatrix>;
  
} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
