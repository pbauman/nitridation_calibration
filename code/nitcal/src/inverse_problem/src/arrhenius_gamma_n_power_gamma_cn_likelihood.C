//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// C++
#include <cmath>

// This class
#include "arrhenius_gamma_n_power_gamma_cn_likelihood.h"

// GRINS
#include "grins/math_constants.h"

// QUESO
#include "uqGslVector.h"
#include "uqGslMatrix.h"

namespace NitridationCalibration
{
  template<class Vec,class Mat>
  ArrheniusGammaNPowerGammaCNLikelihood<Vec,Mat>::ArrheniusGammaNPowerGammaCNLikelihood( int argc,
                                                                                         char** argv,
                                                                                         MPI_Comm mpi_comm,
                                                                                         const GetPot& sip_input,
                                                                                         const GetPot& forward_run_input,
                                                                                         const LikelihoodCommHandler& comm_handler,
                                                                                         const char* prefix, 
                                                                                         const uqVectorSetClass<Vec,Mat>& domain_set,
                                                                                         const bool returns_ln)
  : GammaNGammaCNLikelihoodBase<Vec,Mat>(argc,argv,mpi_comm,
                                         forward_run_input,comm_handler,
                                         prefix,domain_set,returns_ln),
    _gamma0_CN_nom( sip_input( "InverseProblem/gamma_CN_nominal_value", 1.0e-3 ) ),
    _Tref_CN_nom( sip_input( "InverseProblem/Tref_CN_nominal_value", 300.0 ) ),
    _alpha_CN_nom( sip_input( "InverseProblem/alpha_CN_nominal_value", 1.0 ) ),
    _gamma0_N_nom( sip_input( "InverseProblem/gamma0_N_nominal_value", 1.0e-3 ) ),
    _Ta_N_nom( sip_input( "InverseProblem/Ta_N_nominal_value", 300.0 ) )
  {
    return;
  }

  template<class Vec,class Mat>
  ArrheniusGammaNPowerGammaCNLikelihood<Vec,Mat>::~ArrheniusGammaNPowerGammaCNLikelihood()
  {
    return;
  }

  template<class Vec,class Mat>
  void ArrheniusGammaNPowerGammaCNLikelihood<Vec,Mat>::update_parameters( const std::vector<double>& params ) const
  {
    if( params.size() != 5 )
      {
        std::cerr << "Error: params size should be 5 for ArrheniusGammaNPowerGammaCNLikelihood" << std::endl
                  << "       Found size = " << params.size() << std::endl;
        libmesh_error();
      }

    std::vector<libMesh::Real> gamma_CN_params(3);
    gamma_CN_params[0] = params[0]*_gamma0_CN_nom;
    gamma_CN_params[1] = params[1]*_Tref_CN_nom;
    gamma_CN_params[2] = params[2]*_alpha_CN_nom;
    this->_interface.set_gamma_CN_params( gamma_CN_params );

    std::vector<libMesh::Real> gamma_N_params(2);
    gamma_N_params[0] = params[3]*_gamma0_N_nom;
    gamma_N_params[1] = params[4]*_Ta_N_nom;

    this->_interface.set_gamma_N_params( gamma_N_params );

    /*
    if( this->m_env.fullRank() == 0 )
    {
      std::cout << "param = " << params[0] << ", gamma = " << params[0]*_gamma_nom << std::endl;
    }
    */

    return;
  }

  // Instantiate GSL version of this class
  template class ArrheniusGammaNPowerGammaCNLikelihood<uqGslVectorClass,uqGslMatrixClass>;

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
