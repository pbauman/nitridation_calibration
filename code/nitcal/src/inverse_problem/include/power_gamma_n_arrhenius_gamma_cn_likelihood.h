//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#ifndef NITCAL_POWER_GAMMA_N_ARRHENIUS_GAMMA_CN_LIKELIHOOD_H
#define NITCAL_POWER_GAMMA_N_ARRHENIUS_GAMMA_CN_LIKELIHOOD_H

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

#include "gamma_n_gamma_cn_likelihood_base.h"

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  class PowerGammaNArrheniusGammaCNLikelihood : public GammaNGammaCNLikelihoodBase<Vec,Mat>
  {
  public:

    PowerGammaNArrheniusGammaCNLikelihood( int argc,
                                          char** argv,
                                          MPI_Comm mpi_comm,
                                          const GetPot& input,
                                          const char* prefix, 
                                          const uqVectorSetClass<Vec,Mat>& domain_set,
                                          const bool returns_ln );

    virtual ~PowerGammaNArrheniusGammaCNLikelihood();

    virtual void update_parameters( const std::vector<double>& params ) const;

  protected:

    double _gamma0_CN_nom;
    double _Ta_CN_nom;

    double _gamma0_N_nom;
    double _Tref_N_nom;
    double _alpha_N_nom;

  };

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO

#endif // NITCAL_POWER_GAMMA_N_ARRHENIUS_GAMMA_CN_LIKELIHOOD_H
