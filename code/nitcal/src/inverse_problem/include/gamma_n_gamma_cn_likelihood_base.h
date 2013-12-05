//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#ifndef NITCAL_GAMMA_N_GAMMA_CN_LIKELIHOOD_BASE_H
#define NITCAL_GAMMA_N_GAMMA_CN_LIKELIHOOD_BASE_H

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

#include "likelihood_base.h"
#include "mass_loss_helper.h"
#include "average_N_helper.h"

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  class GammaNGammaCNLikelihoodBase : public LikelihoodBase<Vec,Mat>
  {
  public:

    GammaNGammaCNLikelihoodBase( int argc,
                                 char** argv,
                                 MPI_Comm mpi_comm,
                                 const GetPot& forward_run_input,
                                 const LikelihoodCommHandler& comm_handler,
                                 const char* prefix, 
                                 const uqVectorSetClass<Vec,Mat>& domain_set,
                                 const bool returns_ln );

    virtual ~GammaNGammaCNLikelihoodBase();

    virtual double evaluate_likelihood() const;

  protected:

    MassLossHelper _mass_loss;
    
    AverageNHelper _average_n;

  private:

    GammaNGammaCNLikelihoodBase();
    
  };

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO

#endif // NITCAL_GAMMA_N_GAMMA_CN_LIKELIHOOD_BASE_H
