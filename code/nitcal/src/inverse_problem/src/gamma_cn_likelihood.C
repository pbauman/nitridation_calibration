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
#include "gamma_cn_likelihood.h"

// GRINS
#include "grins/math_constants.h"

// QUESO
#include "uqGslVector.h"
#include "uqGslMatrix.h"

namespace NitridationCalibration
{
  template<class Vec,class Mat>
  GammaCNLikelihood<Vec,Mat>::GammaCNLikelihood( int argc,
						 char** argv,
						 MPI_Comm mpi_comm,
						 const std::string& input_filename,
						 const char* prefix, 
						 const uqVectorSetClass<Vec,Mat>& domain_set,
						 const bool returns_ln)
    : LikelihoodBase<Vec,Mat>(prefix,domain_set,returns_ln),
      _interface(argc,argv,mpi_comm, input_filename ),
      _mass_loss(input_filename)
  {
    GetPot input( input_filename );

    _gamma_nom = input( "GammaCNInverseProblem/gamma_nominal_value", 1.0e-3 );

    return;
  }

  template<class Vec,class Mat>
  GammaCNLikelihood<Vec,Mat>::~GammaCNLikelihood()
  {
    return;
  }

  template<class Vec,class Mat>
  void GammaCNLikelihood<Vec,Mat>::update_parameters( const std::vector<double>& params ) const
  {
    libmesh_assert_equal_to( params.size(), 1 );
    
    _interface.set_gamma_CN( params[0]*_gamma_nom );

    if( this->m_env.fullRank() == 0 )
    {
      std::cout << "param = " << params[0] << ", gamma = " << params[0]*_gamma_nom << std::endl;
    }

    return;
  }

  template<class Vec,class Mat>
  double GammaCNLikelihood<Vec,Mat>::evaluate_likelihood() const
  {
    _interface.solve();

    const double computed_mass_loss = _interface.computed_mass_loss();

    // Reset initial guess for next time
    _interface.reset_initial_guess();

    libMesh::Real likelihood_value = _mass_loss.likelihood_value(computed_mass_loss);

    // Now sum over all the data sets distributed across the processors.
    //Parallel::sum(likelihood_value, Parallel::Communicator(this->m_env.subComm().Comm()));

    if( this->m_env.fullRank() == 0 )
    {
      this->_sample_count += 1;

      std::cout << "Sample #" << this->_sample_count << std::endl
      << "mass loss = " << computed_mass_loss << std::endl
      << "log(likelihood) = " << likelihood_value << std::endl << std::endl;
    }

    return likelihood_value;
  }

  // Instantiate GSL version of this class
  template class GammaCNLikelihood<uqGslVectorClass,uqGslMatrixClass>;

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
