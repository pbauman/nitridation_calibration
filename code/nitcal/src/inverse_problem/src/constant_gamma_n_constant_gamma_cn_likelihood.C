//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// C++
#include <cmath>

// This class
#include "constant_gamma_n_constant_gamma_cn_likelihood.h"

// GRINS
#include "grins/math_constants.h"

// QUESO
#include "queso/GslVector.h"
#include "queso/GslMatrix.h"

namespace NitridationCalibration
{
  template<class Vec,class Mat>
  ConstantGammaNConstantGammaCNLikelihood<Vec,Mat>::ConstantGammaNConstantGammaCNLikelihood( int argc,
                                                                                             char** argv,
                                                                                             MPI_Comm mpi_comm,
                                                                                             const GetPot& sip_input,
                                                                                             const GetPot& forward_run_input,
                                                                                             const LikelihoodCommHandler& comm_handler,
                                                                                             const char* prefix, 
                                                                                             const QUESO::VectorSet<Vec,Mat>& domain_set,
                                                                                             const bool returns_ln)
  : GammaNGammaCNLikelihoodBase<Vec,Mat>(argc,argv,mpi_comm,
                                         forward_run_input,comm_handler,
                                         prefix,domain_set,returns_ln),
    _gamma_CN_nom( sip_input( "InverseProblem/gamma_CN_nominal_value", 1.0e-3 ) ),
    _gamma_N_nom( sip_input( "InverseProblem/gamma_N_nominal_value", 1.0e-3 ) )
  {
    return;
  }

  template<class Vec,class Mat>
  ConstantGammaNConstantGammaCNLikelihood<Vec,Mat>::~ConstantGammaNConstantGammaCNLikelihood()
  {
    return;
  }

  template<class Vec,class Mat>
  void ConstantGammaNConstantGammaCNLikelihood<Vec,Mat>::update_parameters( const std::vector<double>& params ) const
  {
    if( params.size() != 2 )
      {
        std::cerr << "Error: params size should be 2 for ConstantGammaNConstantGammaCNLikelihood" << std::endl
                  << "       Found size = " << params.size() << std::endl;
        libmesh_error();
      }
    
    std::vector<libMesh::Real> gamma_CN_params(1);
    gamma_CN_params[0] = params[0]*_gamma_CN_nom;
    this->_interface.set_gamma_CN_params( gamma_CN_params );

    std::vector<libMesh::Real> gamma_N_params(1);
    gamma_N_params[0] = params[1]*_gamma_N_nom;

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
  template class ConstantGammaNConstantGammaCNLikelihood<QUESO::GslVector,QUESO::GslMatrix>;

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
