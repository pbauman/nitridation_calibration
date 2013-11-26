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

// QUESO
#include "uqGslVector.h"
#include "uqGslMatrix.h"

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  LikelihoodBase<Vec,Mat>::LikelihoodBase( int argc, char** argv, MPI_Comm mpi_comm,
                                           const GetPot& input,
                                           const char* prefix, 
                                           const uqVectorSetClass<Vec,Mat>& domain_set,
                                           const bool returns_ln )
    : QuesoLikelihoodInterface<Vec,Mat>( prefix, domain_set, returns_ln ),
      _interface(argc,argv,mpi_comm, input ),
      _sample_count(0)

  {
    return;
  }

  template<class Vec,class Mat>
  LikelihoodBase<Vec,Mat>::~LikelihoodBase()
  {
    return;
  }

  // Instantiate GSL version of this class
  template class LikelihoodBase<uqGslVectorClass,uqGslMatrixClass>;

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
