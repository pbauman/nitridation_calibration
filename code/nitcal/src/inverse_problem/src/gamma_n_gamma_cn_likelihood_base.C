//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// C++
#include <cmath>

// This class
#include "gamma_n_gamma_cn_likelihood_base.h"

// Nitcal
#include "likelihood_comm_handler.h"

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
                                                                     const GetPot& forward_run_input,
                                                                     const LikelihoodCommHandler& comm_handler,
                                                                     const char* prefix, 
                                                                     const uqVectorSetClass<Vec,Mat>& domain_set,
                                                                     const bool returns_ln)
  : LikelihoodBase<Vec,Mat>(argc,argv,mpi_comm,
                            forward_run_input,comm_handler,
                            prefix,domain_set,returns_ln),
    _mass_loss(forward_run_input),
    _average_n(forward_run_input)
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

    /*
    std::cout << std::scientific << std::setprecision(16)
              << "computed_mass_loss = " << computed_mass_loss << std::endl
              << "computed_average_n = " << computed_average_n << std::endl;
    */

    // Reset initial guess for next time
    this->_interface.reset_initial_guess();

    libMesh::Real likelihood_value = 0.0;

    likelihood_value += _mass_loss.likelihood_value(computed_mass_loss);

    likelihood_value += _average_n.likelihood_value(computed_average_n);

    // Now sum over all the data sets distributed across the processors.
    libMesh::Parallel::sum(likelihood_value, libMesh::Parallel::Communicator(this->_comm_handler.get_inter_chain_0_comm()));

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
