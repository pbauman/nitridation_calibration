//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// This class
#include "arrhenius_gamma_n_arrhenius_gamma_cn_sip.h"

// NitCal
#include "arrhenius_gamma_n_arrhenius_gamma_cn_likelihood.h"

// QUESO
#include "uqGslVector.h"
#include "uqGslMatrix.h"

namespace NitridationCalibration
{
  template<class Vec,class Mat>
  ArrheniusGammaNArrheniusGammaCNSIP<Vec,Mat>::ArrheniusGammaNArrheniusGammaCNSIP( uqBaseEnvironmentClass* env,
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
  ArrheniusGammaNArrheniusGammaCNSIP<Vec,Mat>::~ArrheniusGammaNArrheniusGammaCNSIP()
  {
    return;
  }

  template<class Vec,class Mat>
  void ArrheniusGammaNArrheniusGammaCNSIP<Vec,Mat>::create_param_space()
  {
    const unsigned int n_params = 4;

    this->_param_space = new uqVectorSpaceClass<Vec,Mat>( *(this->_queso_env),
							  "param_",
							  n_params,
							  NULL );

    return;
  }

  template<class Vec,class Mat>
  void ArrheniusGammaNArrheniusGammaCNSIP<Vec,Mat>::create_param_domain( const GetPot& input )
  {
    // These are assuming normalized values
    const double gamma0_CN_min = input("InverseProblem/gamma0_CN_min", 0.0);
    const double Ta_CN_min = input("InverseProblem/Ta_CN_min", 0.0);
    const double gamma0_N_min = input("InverseProblem/gamma0_N_min", 0.0);
    const double Ta_N_min = input("InverseProblem/Ta_N_min", 0.0);

    const double gamma0_CN_max = input("InverseProblem/gamma0_CN_max", 1000.0);
    const double Ta_CN_max = input("InverseProblem/Ta_CN_max", 1000.0);
    const double gamma0_N_max = input("InverseProblem/gamma0_N_max", 1000.0);
    const double Ta_N_max = input("InverseProblem/Ta_N_max", 1000.0);

    Vec param_mins( this->_param_space->zeroVector() );
    param_mins[0] = gamma0_CN_min;
    param_mins[1] = Ta_CN_min;
    param_mins[2] = gamma0_N_min;
    param_mins[3] = Ta_N_min;

    Vec param_maxs( this->_param_space->zeroVector() );
    param_maxs[0] = gamma0_CN_max;
    param_maxs[1] = Ta_CN_max;
    param_maxs[2] = gamma0_N_max;
    param_maxs[3] = Ta_N_max;

    this->_param_domain = new uqBoxSubsetClass<Vec,Mat>("param_",
							*(this->_param_space),
							param_mins,
							param_maxs );
    
    return;
  }

  template<class Vec,class Mat>
  void ArrheniusGammaNArrheniusGammaCNSIP<Vec,Mat>::create_likelihood( int argc,
                                                                       char** argv,
                                                                       MPI_Comm mpi_comm,
                                                                       const GetPot& sip_input,
                                                                       const GetPot& forward_run_input )
  {
    this->_likelihood = 
      new ArrheniusGammaNArrheniusGammaCNLikelihood<Vec,Mat>( argc, argv, mpi_comm,
                                                              sip_input,
                                                              forward_run_input,
                                                              *this->_comm_handler.get(),
                                                              "like_",
                                                              *(this->_param_domain),
                                                              true ); // the routine computes [ln(function)]

    return;
  }

  // Instantiate GSL version of this class
  template class ArrheniusGammaNArrheniusGammaCNSIP<uqGslVectorClass,uqGslMatrixClass>;

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
