//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#ifndef NITCAL_CONSTANT_GAMMA_N_CONSTANT_GAMMA_CN_MODEL_H
#define NITCAL_CONSTANT_GAMMA_N_CONSTANT_GAMMA_CN_MODEL_H

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

#include "model_interface_base.h"

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  class ConstantGammaNConstantGammaCNModel : public ModelInterfaceBase<Vec,Mat>
  {
  public:

    ConstantGammaNConstantGammaCNModel( const QUESO::BaseEnvironment& env,
                                        const GetPot& model_input );

    virtual ~ConstantGammaNConstantGammaCNModel(){};

    virtual void update_parameters( const std::vector<double>& param_values,
                                    std::vector<double>& gamma_CN_params,
                                    std::vector<double>& gamma_N_params ) const;

  protected:

    double _gamma_CN_nom;

    double _gamma_N_nom;

  private:

    ConstantGammaNConstantGammaCNModel();

  };

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO

#endif // NITCAL_CONSTANT_GAMMA_N_CONSTANT_GAMMA_CN_MODEL_H
