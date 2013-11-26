//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef NITCAL_GAMMA_CN_BASE_H
#define NITCAL_GAMMA_CN_BASE_H

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

#include "sip_base.h"

namespace NitridationCalibration
{
  template<class Vec,class Mat>
  class ConstantGammaCNSIP : public StatisticalInverseProblemBase<Vec,Mat>
  {
  public:

    ConstantGammaCNSIP( uqBaseEnvironmentClass* env,
                        const std::string& method,
                        int argc,
                        char** argv,
                        const std::string& input_filename );

    virtual ~ConstantGammaCNSIP();

  protected:

    void create_param_space();
    void create_param_domain();
    void create_likelihood( int argc, char** argv, MPI_Comm mpi_comm,
                            const std::string& input_filename );

  };

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO

#endif // NITCAL_GAMMA_CN_BASE_H
