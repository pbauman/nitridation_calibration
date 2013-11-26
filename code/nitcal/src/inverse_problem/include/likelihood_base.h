//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id: likelihood_base.h 38830 2013-04-23 18:54:34Z pbauman $ 
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

#include "queso_likelihood_interface.h"
#include "simulation_interface.h"

#ifndef NITCAL_LIKELIHOOD_BASE_H
#define NITCAL_LIKELIHOOD_BASE_H

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  class LikelihoodBase : public QuesoLikelihoodInterface<Vec,Mat>
  {
  public:

    LikelihoodBase( int argc, char** argv, MPI_Comm mpi_comm,
                    const GetPot& input,
                    const char* prefix, 
                    const uqVectorSetClass<Vec,Mat>& domain_set,
                    const bool returns_ln);

    virtual ~LikelihoodBase();

  protected:

    SimulationInterface _interface;

    // Silly variable to monitor sampling
    // since this problem takes a bit.
    mutable unsigned int _sample_count;

  };

} // end namespace NitridationCalibration

#endif // NITCAL_LIKELIHOOD_BASE_H

#endif // NITCAL_HAVE_QUESO
