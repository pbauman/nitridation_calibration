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
#include "constant_gamma_cn_likelihood.h"

// Nitcal
#include "likelihood_comm_handler.h"

// GRINS
#include "grins/math_constants.h"

// QUESO
#include "queso/GslVector.h"
#include "queso/GslMatrix.h"

namespace NitridationCalibration
{
  template<class Vec,class Mat>
  ConstantGammaCNLikelihood<Vec,Mat>::ConstantGammaCNLikelihood( int argc,
                                                                 char** argv,
                                                                 MPI_Comm mpi_comm,
                                                                 const GetPot& sip_input,
                                                                 const GetPot& forward_run_input,
                                                                 const LikelihoodCommHandler& comm_handler,
                                                                 const char* prefix, 
                                                                 const QUESO::VectorSet<Vec,Mat>& domain_set,
                                                                 const bool returns_ln)
    : LikelihoodBase<Vec,Mat>(argc,argv,mpi_comm,
                              forward_run_input,comm_handler,
                              prefix,domain_set,returns_ln),
    _mass_loss(forward_run_input),
    _gamma_nom( sip_input( "InverseProblem/gamma_CN_nominal_value", 1.0e-3 ) )
  {
    return;
  }

  template<class Vec,class Mat>
  ConstantGammaCNLikelihood<Vec,Mat>::~ConstantGammaCNLikelihood()
  {
    return;
  }

  template<class Vec,class Mat>
  void ConstantGammaCNLikelihood<Vec,Mat>::update_parameters( const std::vector<double>& params ) const
  {
    libmesh_assert_equal_to( params.size(), 1 );
    
    std::vector<libMesh::Real> gamma_CN_params(1);
    gamma_CN_params[0] = params[0]*_gamma_nom;
    this->_interface.set_gamma_CN_params( gamma_CN_params );

    if( this->m_env.fullRank() == 0 )
    {
      std::cout << "param = " << params[0] << ", gamma = " << params[0]*_gamma_nom << std::endl;
    }

    return;
  }

  template<class Vec,class Mat>
  double ConstantGammaCNLikelihood<Vec,Mat>::evaluate_likelihood() const
  {
    this->_interface.solve();

    const double computed_mass_loss = this->_interface.computed_mass_loss();

    // Reset initial guess for next time
    this->_interface.reset_initial_guess();

    libMesh::Real likelihood_value = _mass_loss.likelihood_value(computed_mass_loss);

    // Now sum over all the data sets distributed across the processors.
    if( this->_comm_handler.get_inter0_rank() >= 0 )
      {
        libMesh::Parallel::Communicator comm(this->_comm_handler.get_inter_chain_0_comm());
        comm.sum(likelihood_value);
      }

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
  template class ConstantGammaCNLikelihood<QUESO::GslVector,QUESO::GslMatrix>;

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
