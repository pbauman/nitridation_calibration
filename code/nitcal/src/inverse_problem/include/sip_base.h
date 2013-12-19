//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#ifndef NITCAL_SIP_BASE_H
#define NITCAL_SIP_BASE_H

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// libMesh
#include "libmesh/getpot.h"

// NitCal
#include "queso_sip_interface.h"

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  class StatisticalInverseProblemBase : public QuesoStatisticalInverseProblemInterface<Vec,Mat>
  {
  public:

    StatisticalInverseProblemBase( uqBaseEnvironmentClass* env,
                                   const std::string& method,
                                   int argc,
                                   char** argv,
                                   const std::string& sip_input_filename );

    virtual ~StatisticalInverseProblemBase();

    const GetPot& get_sip_input() const;

    const GetPot& get_forward_run_input() const;

  protected:

    boost::scoped_ptr<GetPot> _sip_input;
    boost::scoped_ptr<GetPot> _forward_run_input;

  private:

    StatisticalInverseProblemBase();

  };

  template<class Vec,class Mat>
  inline
  const GetPot& StatisticalInverseProblemBase<Vec,Mat>::get_sip_input() const
  {
    return (*_sip_input.get());
  }

  template<class Vec,class Mat>
  inline
  const GetPot& StatisticalInverseProblemBase<Vec,Mat>::get_forward_run_input() const
  {
    return (*_forward_run_input.get());
  }

} // end namespace NitridationCalibration
#endif // NITCAL_HAVE_QUESO

#endif // NITCAL_SIP_BASE_H
