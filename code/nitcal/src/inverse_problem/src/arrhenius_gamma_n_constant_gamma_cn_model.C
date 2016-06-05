//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// This class
#include "arrhenius_gamma_n_constant_gamma_cn_model.h"

// QUESO
#include "queso/VectorSet.h"
#include <queso/BoxSubset.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  ArrheniusGammaNConstantGammaCNModel<Vec,Mat>::ArrheniusGammaNConstantGammaCNModel( const QUESO::BaseEnvironment& env,
                                                                                     const GetPot& model_input )
    : ModelInterfaceBase<Vec,Mat>(),
    _Ta_N_nom( model_input( "ModelBounds/Ta_N_nominal_value", 0.0 ) )
  {
    if( !model_input.have_variable("ModelBounds/Ta_N_nominal_value") )
      libmesh_error_msg("ERROR: Could not find input parameter ModelBounds/Ta_N_nominal_value");

    const unsigned int n_params = 3;

    this->_param_space.reset( new QUESO::VectorSpace<Vec,Mat>( env,
                                                               "param_",
                                                               n_params,
                                                               NULL ) );

    if( !model_input.have_variable("ModelBounds/log_gamma_CN_min") )
      libmesh_error_msg("ERROR: Could not find input parameter ModelBounds/log_gamma_CN_min");

    if( !model_input.have_variable("ModelBounds/log_gamma0_N_min") )
      libmesh_error_msg("ERROR: Could not find input parameter ModelBounds/log_gamma0_N_min");

    if( !model_input.have_variable("ModelBounds/Ta_N_min") )
      libmesh_error_msg("ERROR: Could not find input parameter ModelBounds/Ta_N_min");

    if( !model_input.have_variable("ModelBounds/log_gamma_CN_min") )
      libmesh_error_msg("ERROR: Could not find input parameter ModelBounds/log_gamma_CN_max");

    if( !model_input.have_variable("ModelBounds/log_gamma0_N_max") )
      libmesh_error_msg("ERROR: Could not find input parameter ModelBounds/log_gamma0_N_max");

    if( !model_input.have_variable("ModelBounds/Ta_N_min") )
      libmesh_error_msg("ERROR: Could not find input parameter ModelBounds/Ta_N_max");

    // These are assuming normalized values
    const double gamma_CN_min = model_input("ModelBounds/log_gamma_CN_min", 0.0);
    const double gamma0_N_min = model_input("ModelBounds/log_gamma0_N_min", 0.0);
    const double Ta_N_min = model_input("ModelBounds/Ta_N_min", 0.0);

    const double gamma_CN_max = model_input("ModelBounds/log_gamma_CN_max", 1000.0);
    const double gamma0_N_max = model_input("ModelBounds/log_gamma0_N_max", 1000.0);
    const double Ta_N_max = model_input("ModelBounds/Ta_N_max", 1000.0);

    Vec param_mins( this->_param_space->zeroVector() );
    param_mins[0] = gamma_CN_min;
    param_mins[1] = gamma0_N_min;
    param_mins[2] = Ta_N_min;

    Vec param_maxs( this->_param_space->zeroVector() );
    param_maxs[0] = gamma_CN_max;
    param_maxs[1] = gamma0_N_max;
    param_maxs[2] = Ta_N_max;

    this->_param_domain.reset( new QUESO::BoxSubset<Vec,Mat>( "param_",
                                                              *(this->_param_space),
                                                              param_mins,
                                                              param_maxs ) );
  }

  template<class Vec,class Mat>
  void ArrheniusGammaNConstantGammaCNModel<Vec,Mat>::update_parameters( const std::vector<double>& param_values,
                                                                        std::vector<double>& gamma_CN_params,
                                                                        std::vector<double>& gamma_N_params ) const
  {
    if( param_values.size() != 3 )
      {
        std::cerr << "Error: params size should be 3 for ArrheniusGammaNConstantGammaCNModel" << std::endl
                  << "       Found size = " << param_values.size() << std::endl;
        libmesh_error();
      }

    gamma_CN_params.resize(1);
    gamma_CN_params[0] = std::pow(10, param_values[0]);

    gamma_N_params.resize(2);
    gamma_N_params[0] = std::pow(10, param_values[1]);
    gamma_N_params[1] = param_values[2]*_Ta_N_nom;
  }

  // Instantiate GSL version of this class
  template class ArrheniusGammaNConstantGammaCNModel<QUESO::GslVector,QUESO::GslMatrix>;

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
