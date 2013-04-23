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
#include "sip_base.h"

// QUESO
#include "uqGslVector.h"
#include "uqGslMatrix.h"

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  StatisticalInverseProblemBase<Vec,Mat>::StatisticalInverseProblemBase( uqBaseEnvironmentClass *env,
									 const std::string method )
    : _queso_env(env),
      _method(method)
  {
    return;
  }

  template<class Vec,class Mat>
  StatisticalInverseProblemBase<Vec,Mat>::~StatisticalInverseProblemBase( )
  {
    delete _ip;
    
    return;
  }

  template<class Vec,class Mat>
  void StatisticalInverseProblemBase<Vec,Mat>::create_sip( )
  {
    _ip = new uqStatisticalInverseProblemClass<Vec,Mat>("", // No extra prefix before the default "ip_" prefix
							NULL,
							*_prior,
							*_likelihood,
							*_posterior);
    
    return;
  }
  
  template<class Vec,class Mat>
  void StatisticalInverseProblemBase<Vec,Mat>::solve( )
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
    const uqBaseVectorRealizerClass<Vec,Mat>& postRealizer = postRv.realizer();
    uqVectorSpaceClass<Vec,Mat> auxSpace(env,"",postRealizer.subPeriod(),NULL);
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

  /* ------------------------- Instantiate -------------------------*/
  template class StatisticalInverseProblemBase<uqGslVectorClass,uqGslMatrixClass>;
  
} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
