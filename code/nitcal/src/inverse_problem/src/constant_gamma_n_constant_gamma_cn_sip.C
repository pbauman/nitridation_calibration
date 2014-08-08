//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// This class
#include "constant_gamma_n_constant_gamma_cn_sip.h"

// NitCal
#include "constant_gamma_n_constant_gamma_cn_likelihood.h"

// QUESO
#include "queso/GslVector.h"
#include "queso/GslMatrix.h"

// libMesh
#include "libmesh/getpot.h"

namespace NitridationCalibration
{
  template<class Vec,class Mat>
  ConstantGammaNConstantGammaCNSIP<Vec,Mat>::ConstantGammaNConstantGammaCNSIP( QUESO::BaseEnvironment* env,
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

    return;
  }

  template<class Vec,class Mat>
  ConstantGammaNConstantGammaCNSIP<Vec,Mat>::~ConstantGammaNConstantGammaCNSIP()
  {
    return;
  }

  template<class Vec,class Mat>
  void ConstantGammaNConstantGammaCNSIP<Vec,Mat>::create_param_space()
  {
    const unsigned int n_params = 2;

    this->_param_space = new QUESO::VectorSpace<Vec,Mat>( *(this->_queso_env),
							  "param_",
							  n_params,
							  NULL );

    return;
  }

  template<class Vec,class Mat>
  void ConstantGammaNConstantGammaCNSIP<Vec,Mat>::create_param_domain( const GetPot& input )
  {
    // These are assuming normalized values
    const double gamma_CN_min = input("InverseProblem/gamma_CN_min", 0.0);
    const double gamma_N_min = input("InverseProblem/gamma_N_min", 0.0);

    const double gamma_CN_max = input("InverseProblem/gamma_CN_max", 1000.0);
    const double gamma_N_max = input("InverseProblem/gamma_N_max", 1000.0);

    Vec param_mins( this->_param_space->zeroVector() );
    param_mins[0] = gamma_CN_min;
    param_mins[1] = gamma_N_min;

    Vec param_maxs( this->_param_space->zeroVector() );
    param_maxs[0] = gamma_CN_max;
    param_maxs[1] = gamma_N_max;

    this->_param_domain = new QUESO::BoxSubset<Vec,Mat>("param_",
							*(this->_param_space),
							param_mins,
							param_maxs );
    
    return;
  }

  template<class Vec,class Mat>
  void ConstantGammaNConstantGammaCNSIP<Vec,Mat>::create_likelihood( int argc,
                                                                     char** argv,
                                                                     MPI_Comm mpi_comm,
                                                                     const GetPot& sip_input,
                                                                     const GetPot& forward_run_input )
  {
    this->_likelihood = 
      new ConstantGammaNConstantGammaCNLikelihood<Vec,Mat>( argc, argv, mpi_comm,
                                                            sip_input,
                                                            forward_run_input,
                                                            *this->_comm_handler.get(),
                                                            "like_",
                                                            *(this->_param_domain),
                                                            true ); // the routine computes [ln(function)]

    return;
  }

  // Instantiate GSL version of this class
  template class ConstantGammaNConstantGammaCNSIP<QUESO::GslVector,QUESO::GslMatrix>;

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
