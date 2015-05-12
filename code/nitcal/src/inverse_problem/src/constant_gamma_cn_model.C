//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// This class
#include "constant_gamma_cn_model.h"

// QUESO
#include "queso/VectorSet.h"
#include <queso/BoxSubset.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  ConstantGammaCNModel<Vec,Mat>::ConstantGammaCNModel( int argc, char** argv,
                                                       const QUESO::BaseEnvironment& queso_env,
                                                       const GetPot& forward_run_input,
                                                       const GetPot& model_input )
    : ModelInterfaceBase<Vec,Mat>(argc,argv,queso_env,forward_run_input),
    _gamma_nom(model_input( "ModelBounds/gamma_CN_nominal_value", 1.0e-3 ))
  {
    const unsigned int n_params = 1;

    this->_param_space.reset( new QUESO::VectorSpace<Vec,Mat>( this->_queso_env,
                                                               "param_",
                                                               n_params,
                                                               NULL ) );

    // These are assuming normalized values
    const double min = model_input("ModelBounds/gamma_CN_min", 0.0);
    const double max = model_input("ModelBounds/gamma_CN_max", 1000.0);

    Vec param_mins( this->_param_space->zeroVector() );
    param_mins[0] = min;

    Vec param_maxs( this->_param_space->zeroVector() );
    param_maxs[0] = max;

    this->_param_domain.reset( new QUESO::BoxSubset<Vec,Mat>( "param_",
                                                              *(this->_param_space),
                                                              param_mins,
                                                              param_maxs ) );
  }

  template<class Vec,class Mat>
  void ConstantGammaCNModel<Vec,Mat>::update_parameters( const std::vector<double>& param_values )
  {
    libmesh_assert_equal_to( param_values.size(), 1 );

    std::vector<libMesh::Real> gamma_CN_params(1);

    gamma_CN_params[0] = param_values[0]*_gamma_nom;

    this->_interface.set_gamma_CN_params( gamma_CN_params );

    if( this->_queso_env.fullRank() == 0 )
    {
      std::cout << "param = " << param_values[0]
                << ", gamma = " << param_values[0]*_gamma_nom << std::endl;
    }
  }

  // Instantiate GSL version of this class
  template class ConstantGammaCNModel<QUESO::GslVector,QUESO::GslMatrix>;

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
