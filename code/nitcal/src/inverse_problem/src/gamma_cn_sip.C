//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// This class
#include "gamma_cn_sip.h"

// NitCal
#include "gamma_cn_likelihood.h"

// QUESO
#include "uqGslVector.h"
#include "uqGslMatrix.h"

namespace NitridationCalibration
{
  template<class Vec,class Mat>
  GammaCNSIP<Vec,Mat>::GammaCNSIP( uqBaseEnvironmentClass* env,
				   const std::string& method,
				   int argc,
				   char** argv,
				   const std::string& libMesh_input_filename )
    : StatisticalInverseProblemBase<Vec,Mat>( env, method )
  {
    this->create_param_space();
    this->create_param_domain();
    this->create_prior();
    this->create_posterior();
    this->create_likelihood(argc,argv,env->subComm().Comm(),libMesh_input_filename);

    this->create_sip();
    
    return;
  }

  template<class Vec,class Mat>
  GammaCNSIP<Vec,Mat>::~GammaCNSIP()
  {
    return;
  }

  template<class Vec,class Mat>
  void GammaCNSIP<Vec,Mat>::create_param_space()
  {
    const unsigned int n_params = 1;

    this->_param_space = new uqVectorSpaceClass<Vec,Mat>( *(this->_queso_env),
							  "param_",
							  n_params,
							  NULL );

    return;
  }

  template<class Vec,class Mat>
  void GammaCNSIP<Vec,Mat>::create_param_domain()
  {
    // These are assuming normalized values
    const double min = 0.0;
    const double max = 10.0;

    Vec param_mins( this->_param_space->zeroVector() );
    param_mins[0] = min;

    Vec param_maxs( this->_param_space->zeroVector() );
    param_maxs[0] = max;

    this->_param_domain = new uqBoxSubsetClass<Vec,Mat>("param_",
							*(this->_param_space),
							param_mins,
							param_maxs );
    
    return;
  }

  template<class Vec,class Mat>
  void GammaCNSIP<Vec,Mat>::create_prior()
  {
    this->_prior =
      new uqUniformVectorRVClass<Vec,Mat>("prior_", // Extra prefix before the default "rv_"
					  *(this->_param_domain) );

    return;
  }

  template<class Vec,class Mat>
  void GammaCNSIP<Vec,Mat>::create_posterior()
  {
    this->_posterior =
      new uqGenericVectorRVClass<Vec,Mat>("post_", // Extra prefix before the default "rv_" prefix
					  *(this->_param_space) );

    return;
  }

  template<class Vec,class Mat>
  void GammaCNSIP<Vec,Mat>::create_likelihood(int argc,
					      char** argv,
					      MPI_Comm mpi_comm,
					      const std::string& input_filename)
  {
    this->_likelihood = 
      new GammaCNLikelihood<Vec,Mat>(argc, argv, mpi_comm, input_filename,
				     "like_",
				     *(this->_param_domain),
				     true ); // the routine computes [ln(function)]

    return;
  }

  // Instantiate GSL version of this class
  template class GammaCNSIP<uqGslVectorClass,uqGslMatrixClass>;

} // end namespace NitridationCalibration
