//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// This class
#include "arrhenius_gamma_n_arrhenius_gamma_cn_model.h"

// QUESO
#include "queso/VectorSet.h"
#include <queso/BoxSubset.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  ArrheniusGammaNArrheniusGammaCNModel<Vec,Mat>::ArrheniusGammaNArrheniusGammaCNModel( const QUESO::BaseEnvironment& env,
                                                                                       const GetPot& model_input )
  : ModelInterfaceBase<Vec,Mat>(),
    _gamma0_CN_nom( model_input( "ModelBounds/gamma0_CN_nominal_value", 0.0 ) ),
    _Ta_CN_nom( model_input( "ModelBounds/Ta_CN_nominal_value", 0.0 ) ),
    _gamma0_N_nom( model_input( "ModelBounds/gamma0_N_nominal_value", 0.0 ) ),
    _Ta_N_nom( model_input( "ModelBounds/Ta_N_nominal_value", 0.0 ) )
  {
    if( !input.have_variable("ModelBounds/gamma0_CN_nominal_value") )
      libmesh_error_msg("ERROR: Could not find input parameter "+"ModelBounds/gamma0_CN_nominal_value");

    if( !input.have_variable("ModelBounds/Ta_CN_nominal_value") )
      libmesh_error_msg("ERROR: Could not find input parameter "+"ModelBounds/Ta_CN_nominal_value");

    if( !input.have_variable("ModelBounds/gamma0_N_nominal_value") )
      libmesh_error_msg("ERROR: Could not find input parameter "+"ModelBounds/gamma0_N_nominal_value");

    if( !input.have_variable("ModelBounds/Ta_N_nominal_value") )
      libmesh_error_msg("ERROR: Could not find input parameter "+"ModelBounds/Ta_N_nominal_value");

    const unsigned int n_params = 4;

    this->_param_space.reset( new QUESO::VectorSpace<Vec,Mat>( env,
                                                               "param_",
                                                               n_params,
                                                               NULL ) );

    const double gamma0_CN_min = model_input("ModelBounds/gamma0_CN_min", 0.0);
    const double Ta_CN_min = model_input("ModelBounds/Ta_CN_min", 0.0);
    const double gamma0_N_min = model_input("ModelBounds/gamma0_N_min", 0.0);
    const double Ta_N_min = model_input("ModelBounds/Ta_N_min", 0.0);

    const double gamma0_CN_max = model_input("ModelBounds/gamma0_CN_max", 1000.0);
    const double Ta_CN_max = model_input("ModelBounds/Ta_CN_max", 1000.0);
    const double gamma0_N_max = model_input("ModelBounds/gamma0_N_max", 1000.0);
    const double Ta_N_max = model_input("ModelBounds/Ta_N_max", 1000.0);

    Vec param_mins( this->_param_space->zeroVector() );
    param_mins[0] = gamma0_CN_min;
    param_mins[1] = Ta_CN_min;
    param_mins[2] = gamma0_N_min;
    param_mins[3] = Ta_N_min;

    Vec param_maxs( this->_param_space->zeroVector() );
    param_maxs[0] = gamma0_CN_max;
    param_maxs[1] = Ta_CN_max;
    param_maxs[2] = gamma0_N_max;
    param_maxs[3] = Ta_N_max;

    this->_param_domain.reset( new QUESO::BoxSubset<Vec,Mat>( "param_",
                                                              *(this->_param_space),
                                                              param_mins,
                                                              param_maxs ) );
  }

  template<class Vec,class Mat>
  void ArrheniusGammaNArrheniusGammaCNModel<Vec,Mat>::update_parameters( const std::vector<double>& param_values,
                                                                         std::vector<double>& gamma_CN_params,
                                                                         std::vector<double>& gamma_N_params ) const
  {
    if( param_values.size() != 4 )
      {
        std::cerr << "Error: params size should be 4 for ArrheniusGammaNArrheniusGammaCNModel" << std::endl
                  << "       Found size = " << param_values.size() << std::endl;
        libmesh_error();
      }

    gamma_CN_params.resize(2);
    gamma_CN_params[0] = param_values[0]*_gamma0_CN_nom;
    gamma_CN_params[1] = param_values[1]*_Ta_CN_nom;

    gamma_N_params.resize(2);
    gamma_N_params[0] = param_values[2]*_gamma0_N_nom;
    gamma_N_params[1] = param_values[3]*_Ta_N_nom;
  }

  // Instantiate GSL version of this class
  template class ArrheniusGammaNArrheniusGammaCNModel<QUESO::GslVector,QUESO::GslMatrix>;

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
