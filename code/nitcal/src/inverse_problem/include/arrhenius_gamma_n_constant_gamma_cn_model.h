//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#ifndef NITCAL_ARRHENIUS_GAMMA_N_CONSTANT_GAMMA_CN_MODEL_H
#define NITCAL_ARRHENIUS_GAMMA_N_CONSTANT_GAMMA_CN_MODEL_H

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

#include "model_interface_base.h"

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  class ArrheniusGammaNConstantGammaCNModel : public ModelInterfaceBase<Vec,Mat>
  {
  public:

    ArrheniusGammaNConstantGammaCNModel( int argc, char** argv,
                                        const QUESO::BaseEnvironment& queso_env,
                                        const GetPot& forward_run_input,
                                        const GetPot& model_input );

    virtual ~ArrheniusGammaNConstantGammaCNModel(){};

  protected:

    virtual void update_parameters( const std::vector<double>& param_values );

    double _gamma_CN_nom;

    double _gamma0_N_nom;
    double _Ta_N_nom;

  private:

    ArrheniusGammaNConstantGammaCNModel();

  };

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO

#endif // NITCAL_ARRHENIUS_GAMMA_N_CONSTANT_GAMMA_CN_MODEL_H
