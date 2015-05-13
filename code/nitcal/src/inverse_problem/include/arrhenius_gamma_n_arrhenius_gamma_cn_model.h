//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#ifndef NITCAL_ARRHENIUS_GAMMA_N_ARRHENIUS_GAMMA_CN_MODEL_H
#define NITCAL_ARRHENIUS_GAMMA_N_ARRHENIUS_GAMMA_CN_MODEL_H

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

#include "model_interface_base.h"

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  class ArrheniusGammaNArrheniusGammaCNModel : public ModelInterfaceBase<Vec,Mat>
  {
  public:

    ArrheniusGammaNArrheniusGammaCNModel( const QUESO::BaseEnvironment& env,
                                          const GetPot& model_input );

    virtual ~ArrheniusGammaNArrheniusGammaCNModel(){};

    virtual void update_parameters( const std::vector<double>& param_values,
                                    std::vector<double>& gamma_CN_params,
                                    std::vector<double>& gamma_N_params ) const;

  protected:

    double _gamma0_CN_nom;
    double _Ta_CN_nom;

    double _gamma0_N_nom;
    double _Ta_N_nom;

  private:

    ArrheniusGammaNArrheniusGammaCNModel();

  };

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO

#endif // NITCAL_ARRHENIUS_GAMMA_N_ARRHENIUS_GAMMA_CN_MODEL_H
