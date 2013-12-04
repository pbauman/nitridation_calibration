//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// This class
#include "constant_gamma_cn_sip.h"

// NitCal
#include "constant_gamma_cn_likelihood.h"

// QUESO
#include "uqGslVector.h"
#include "uqGslMatrix.h"

// libMesh
#include "libmesh/getpot.h"

namespace NitridationCalibration
{
  template<class Vec,class Mat>
  ConstantGammaCNSIP<Vec,Mat>::ConstantGammaCNSIP( uqBaseEnvironmentClass* env,
                                                   const std::string& method,
                                                   int argc,
                                                   char** argv,
                                                   const std::string& sip_input_filename )
    : StatisticalInverseProblemBase<Vec,Mat>(env, method, argc, argv, sip_input_filename)
  {
    this->create_param_space();
    this->create_param_domain( (*(this->_sip_input.get())) );
    this->create_prior();
    this->create_posterior();
    this->create_likelihood( argc, argv, this->_comm_handler->get_split_chain_comm(),
                             (*(this->_sip_input.get())),
                             (*(this->_forward_run_input.get())) );

    this->create_sip();
    
    // We don't need these anymore
    this->_sip_input.reset();
    this->_forward_run_input.reset();

    return;
  }

  template<class Vec,class Mat>
  ConstantGammaCNSIP<Vec,Mat>::~ConstantGammaCNSIP()
  {
    return;
  }

  template<class Vec,class Mat>
  void ConstantGammaCNSIP<Vec,Mat>::create_param_space()
  {
    const unsigned int n_params = 1;

    this->_param_space = new uqVectorSpaceClass<Vec,Mat>( *(this->_queso_env),
							  "param_",
							  n_params,
							  NULL );

    return;
  }

  template<class Vec,class Mat>
  void ConstantGammaCNSIP<Vec,Mat>::create_param_domain(const GetPot& input)
  {
    // These are assuming normalized values
    const double min = input("InverseProblem/gamma_CN_min", 0.0);
    const double max = input("InverseProblem/gamma_CN_max", 1000.0);

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
  void ConstantGammaCNSIP<Vec,Mat>::create_likelihood( int argc, char** argv, MPI_Comm mpi_comm,
                                                       const GetPot& sip_input,
                                                       const GetPot& forward_run_input )
  {
    this->_likelihood = 
      new ConstantGammaCNLikelihood<Vec,Mat>(argc, argv, mpi_comm,
                                             sip_input,
                                             forward_run_input,
                                             *this->_comm_handler.get(),
                                             "like_",
                                             *(this->_param_domain),
                                             true ); // the routine computes [ln(function)]

    return;
  }

  // Instantiate GSL version of this class
  template class ConstantGammaCNSIP<uqGslVectorClass,uqGslMatrixClass>;

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
