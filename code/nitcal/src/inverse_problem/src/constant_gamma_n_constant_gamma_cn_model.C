//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// This class
#include "constant_gamma_n_constant_gamma_cn_model.h"

// QUESO
#include "queso/VectorSet.h"
#include <queso/BoxSubset.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  ConstantGammaNConstantGammaCNModel<Vec,Mat>::ConstantGammaNConstantGammaCNModel( int argc, char** argv,
                                                                                   const QUESO::BaseEnvironment& queso_env,
                                                                                   const GetPot& forward_run_input,
                                                                                   const GetPot& model_input )
    : ModelInterfaceBase<Vec,Mat>(argc,argv,queso_env,forward_run_input),
    _gamma_CN_nom( model_input( "ModelBounds/gamma_CN_nominal_value", 1.0e-3 ) ),
    _gamma_N_nom( model_input( "ModelBounds/gamma_N_nominal_value", 1.0e-3 ) )
  {
    const unsigned int n_params = 2;

    this->_param_space.reset( new QUESO::VectorSpace<Vec,Mat>( this->_queso_env,
                                                               "param_",
                                                               n_params,
                                                               NULL ) );

    // These are assuming normalized values
    const double gamma_CN_min = model_input("ModelBounds/gamma_CN_min", 0.0);
    const double gamma_N_min = model_input("ModelBounds/gamma_N_min", 0.0);

    const double gamma_CN_max = model_input("ModelBounds/gamma_CN_max", 1000.0);
    const double gamma_N_max = model_input("ModelBounds/gamma_N_max", 1000.0);

    Vec param_mins( this->_param_space->zeroVector() );
    param_mins[0] = gamma_CN_min;
    param_mins[1] = gamma_N_min;

    Vec param_maxs( this->_param_space->zeroVector() );
    param_maxs[0] = gamma_CN_max;
    param_maxs[1] = gamma_N_max;

    this->_param_domain.reset( new QUESO::BoxSubset<Vec,Mat>( "param_",
                                                              *(this->_param_space),
                                                              param_mins,
                                                              param_maxs ) );
  }

  template<class Vec,class Mat>
  void ConstantGammaNConstantGammaCNModel<Vec,Mat>::update_parameters( const std::vector<double>& param_values )
  {
    if( param_values.size() != 2 )
      {
        std::cerr << "Error: params size should be 2 for ConstantGammaNConstantGammaCNModel" << std::endl
                  << "       Found size = " << param_values.size() << std::endl;
        libmesh_error();
      }

    std::vector<libMesh::Real> gamma_CN_params(1);
    gamma_CN_params[0] = param_values[0]*_gamma_CN_nom;
    this->_interface.set_gamma_CN_params( gamma_CN_params );

    std::vector<libMesh::Real> gamma_N_params(1);
    gamma_N_params[0] = param_values[1]*_gamma_N_nom;

    this->_interface.set_gamma_N_params( gamma_N_params );
  }

  // Instantiate GSL version of this class
  template class ConstantGammaNConstantGammaCNModel<QUESO::GslVector,QUESO::GslMatrix>;

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
