//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id: likelihood_base.h 38830 2013-04-23 18:54:34Z pbauman $ 
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// This class
#include "likelihood_base.h"

// NitCal
#include "likelihood_comm_handler.h"

// QUESO
#include "queso/GslVector.h"
#include "queso/GslMatrix.h"

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  LikelihoodBase<Vec,Mat>::LikelihoodBase( int argc, char** argv, MPI_Comm mpi_comm,
                                           const GetPot& forward_run_input,
                                           const LikelihoodCommHandler& comm_handler,
                                           const char* prefix, 
                                           const QUESO::VectorSet<Vec,Mat>& domain_set,
                                           const bool returns_ln )
    : QuesoLikelihoodInterface<Vec,Mat>( prefix, domain_set, returns_ln ),
    _interface(argc,argv,mpi_comm, forward_run_input ),
    _sample_count(0),
    _comm_handler(comm_handler)

  {
    return;
  }

  template<class Vec,class Mat>
  LikelihoodBase<Vec,Mat>::~LikelihoodBase()
  {
    return;
  }

  // Instantiate GSL version of this class
  template class LikelihoodBase<QUESO::GslVector,QUESO::GslMatrix>;

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
