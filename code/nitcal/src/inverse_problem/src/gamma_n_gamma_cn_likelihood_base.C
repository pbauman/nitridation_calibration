//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// C++
#include <cmath>

// This class
#include "gamma_n_gamma_cn_likelihood_base.h"

// GRINS
#include "grins/math_constants.h"

// QUESO
#include "uqGslVector.h"
#include "uqGslMatrix.h"

namespace NitridationCalibration
{
  template<class Vec,class Mat>
  GammaNGammaCNLikelihoodBase<Vec,Mat>::GammaNGammaCNLikelihoodBase( int argc,
                                                                     char** argv,
                                                                     MPI_Comm mpi_comm,
                                                                     const std::string& input_filename,
                                                                     const char* prefix, 
                                                                     const uqVectorSetClass<Vec,Mat>& domain_set,
                                                                     const bool returns_ln)
  : LikelihoodBase<Vec,Mat>(argc,argv,mpi_comm,input_filename,prefix,domain_set,returns_ln),
    _mass_loss(input_filename),
    _average_n(input_filename)
  {
    return;
  }

  template<class Vec,class Mat>
  GammaNGammaCNLikelihoodBase<Vec,Mat>::~GammaNGammaCNLikelihoodBase()
  {
    return;
  }

  template<class Vec,class Mat>
  double GammaNGammaCNLikelihoodBase<Vec,Mat>::evaluate_likelihood() const
  {
    this->_interface.solve();

    const double computed_mass_loss = this->_interface.computed_mass_loss();

    const double computed_average_n = this->_interface.computed_average_n();

    // Reset initial guess for next time
    this->_interface.reset_initial_guess();

    libMesh::Real likelihood_value = 0.0;

    likelihood_value += _mass_loss.likelihood_value(computed_mass_loss);

    likelihood_value += _average_n.likelihood_value(computed_average_n);

    // Now sum over all the data sets distributed across the processors.
    //Parallel::sum(likelihood_value, Parallel::Communicator(this->m_env.subComm().Comm()));

    if( this->m_env.fullRank() == 0 )
      {
        this->_sample_count += 1;

        std::cout << "Sample #" << this->_sample_count << std::endl;
      }

    return likelihood_value;
  }

  // Instantiate GSL version of this class
  template class GammaNGammaCNLikelihoodBase<uqGslVectorClass,uqGslMatrixClass>;

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
