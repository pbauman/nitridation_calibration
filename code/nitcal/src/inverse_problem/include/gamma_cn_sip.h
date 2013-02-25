//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef NITCAL_GAMMA_CN_BASE_H
#define NITCAL_GAMMA_CN_BASE_H

#include "sip_base.h"

namespace NitridationCalibration
{
  template<class Vec,class Mat>
  class GammaCNSIP : public StatisticalInverseProblemBase<Vec,Mat>
  {
  public:

    GammaCNSIP( uqBaseEnvironmentClass* env,
		const std::string& method,
		int argc,
		char** argv,
		const std::string& input_filename );

    virtual ~GammaCNSIP();

  protected:

    virtual void create_param_space();
    virtual void create_param_domain();
    virtual void create_prior();
    virtual void create_posterior();
    virtual void create_likelihood(int argc,
				   char** argv,
				   MPI_Comm mpi_comm,
				   const std::string& input_filename);

  };

} // end namespace NitridationCalibration

#endif // NITCAL_GAMMA_CN_BASE_H
