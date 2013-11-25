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
#include "uqGslVector.h"
#include "uqGslMatrix.h"

namespace NitridationCalibration
{
  template<class Vec,class Mat>
  ConstantGammaNConstantGammaCNLikelihood<Vec,Mat>::ConstantGammaNConstantGammaCNLikelihood( int argc,
                                                                                             char** argv,
                                                                                             MPI_Comm mpi_comm,
                                                                                             const std::string& input_filename,
                                                                                             const char* prefix, 
                                                                                             const uqVectorSetClass<Vec,Mat>& domain_set,
                                                                                             const bool returns_ln)
    : GammaNGammaCNLikelihoodBase<Vec,Mat>(argc,argv,mpi_comm,input_filename,prefix,domain_set,returns_ln)
  {
    GetPot input( input_filename );

    _gamma_CN_nom = input( "InverseProblem/gamma_CN_nominal_value", 1.0e-3 );

    _gamma_N_nom = input( "InverseProblem/gamma_N_nominal_value", 1.0e-3 );

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
    
    this->_interface.set_gamma_CN( params[0]*_gamma_CN_nom );

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
  template class ConstantGammaNConstantGammaCNLikelihood<uqGslVectorClass,uqGslMatrixClass>;

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
