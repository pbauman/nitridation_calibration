//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// This class
#include "model_builder.h"

// NitCal
#include <constant_gamma_n_constant_gamma_cn_model.h>
#include <arrhenius_gamma_n_constant_gamma_cn_model.h>
#include <arrhenius_gamma_n_arrhenius_gamma_cn_model.h>

namespace NitridationCalibration
{
  template<class Vec,class Mat>
  ModelInterfaceBase<Vec,Mat>* ModelBuilder<Vec,Mat>::build_model( const QUESO::BaseEnvironment& queso_env,
                                                                   const GetPot& model_input )
  {
    ModelInterfaceBase<Vec,Mat>* model_ptr = NULL;

    std::string model_type = model_input( "ModelType/model", "DIE!" );

    if( model_type == std::string("constant_gamma_n_constant_gamma_cn") )
      model_ptr = new NitridationCalibration::ConstantGammaNConstantGammaCNModel<QUESO::GslVector,QUESO::GslMatrix>(queso_env, model_input);

    else if( model_type == std::string("arrhenius_gamma_n_constant_gamma_cn") )
      model_ptr = new NitridationCalibration::ArrheniusGammaNConstantGammaCNModel<QUESO::GslVector,QUESO::GslMatrix>(queso_env, model_input);

    else if( model_type == std::string("arrhenius_gamma_n_constant_gamma_cn") )
      model_ptr = new NitridationCalibration::ArrheniusGammaNArrheniusGammaCNModel<QUESO::GslVector,QUESO::GslMatrix>(queso_env, model_input);

    else
      {
        std::cerr << "Error: Invalid Model type! Found " << model_type << std::endl;
        MPI_Abort(MPI_COMM_WORLD,1);
      }

    return model_ptr;
  }

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
