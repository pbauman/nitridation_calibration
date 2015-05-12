//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#ifndef NITCAL_CONSTANT_GAMMA_CN_MODEL_H
#define NITCAL_CONSTANT_GAMMA_CN_MODEL_H

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

#include "model_interface_base.h"

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  class ConstantGammaCNModel : public ModelInterfaceBase<Vec,Mat>
  {
  public:

    ConstantGammaCNModel( int argc, char** argv,
                          const QUESO::BaseEnvironment& queso_env,
                          const GetPot& forward_run_input,
                          const GetPot& model_input );

    virtual ~ConstantGammaCNModel(){};

  protected:

    virtual void update_parameters( const std::vector<double>& param_values );

    double _gamma_nom;

  };

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO

#endif // NITCAL_CONSTANT_GAMMA_CN_MODEL_H
