//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#ifndef NITCAL_ARRHENIUS_GAMMA_N_POWER_GAMMA_CN_SIP_H
#define NITCAL_ARRHENIUS_GAMMA_N_POWER_GAMMA_CN_SIP_H

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

#include "sip_base.h"

class GetPot;

namespace NitridationCalibration
{
  template<class Vec,class Mat>
  class ArrheniusGammaNPowerGammaCNSIP : public StatisticalInverseProblemBase<Vec,Mat>
  {
  public:

    ArrheniusGammaNPowerGammaCNSIP( uqBaseEnvironmentClass* env,
                                      const std::string& method,
                                      int argc,
                                      char** argv,
                                      const std::string& input_filename );

    virtual ~ArrheniusGammaNPowerGammaCNSIP();

  protected:

    void create_param_space();
    void create_param_domain( const GetPot& input );
    void create_likelihood( int argc, char** argv, MPI_Comm mpi_comm,
                            const GetPot& sip_input,
                            const GetPot& forward_run_input );

  };

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO

#endif // NITCAL_ARRHENIUS_GAMMA_N_POWER_GAMMA_CN_SIP_H