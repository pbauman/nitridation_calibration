//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#ifndef NITCAL_CONSTANT_GAMMA_N_CONSTANT_GAMMA_CN_LIKELIHOOD_H
#define NITCAL_CONSTANT_GAMMA_N_CONSTANT_GAMMA_CN_LIKELIHOOD_H

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

#include "gamma_n_gamma_cn_likelihood_base.h"

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  class ConstantGammaNConstantGammaCNLikelihood : public GammaNGammaCNLikelihoodBase<Vec,Mat>
  {
  public:

    ConstantGammaNConstantGammaCNLikelihood( int argc,
                                             char** argv,
                                             MPI_Comm mpi_comm,
                                             const GetPot& sip_input,
                                             const GetPot& forward_run_input,
                                             const LikelihoodCommHandler& comm_handler,
                                             const char* prefix, 
                                             const QUESO::VectorSet<Vec,Mat>& domain_set,
                                             const bool returns_ln );

    virtual ~ConstantGammaNConstantGammaCNLikelihood();

    virtual void update_parameters( const std::vector<double>& params ) const;

  protected:

    double _gamma_CN_nom;

    double _gamma_N_nom;

  private:

    ConstantGammaNConstantGammaCNLikelihood();

  };

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO

#endif // NITCAL_CONSTANT_GAMMA_N_CONSTANT_GAMMA_CN_LIKELIHOOD_H
