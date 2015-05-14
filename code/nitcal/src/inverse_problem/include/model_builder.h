//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

#include "model_interface_base.h"

namespace NitridationCalibration
{
  template<class Vec,class Mat>
  class ModelBuilder
  {
  public:
    static ModelInterfaceBase<Vec,Mat>* build_model( const QUESO::BaseEnvironment& queso_env,
                                                     const GetPot& model_input );
  };
} // end namespace NitridationCalibration


#endif // NITCAL_HAVE_QUESO
