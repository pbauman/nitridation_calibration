//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// This class
#include "full_model_evaluator.h"

// QUESO
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  FullModelEvaluator<Vec,Mat>::FullModelEvaluator( int argc, char** argv,
                                                   const QUESO::BaseEnvironment& queso_env,
                                                   const GetPot& forward_run_input,
                                                   MPI_Comm model_comm,
                                                   const ModelInterfaceBase<Vec,Mat>& model )
    : ModelEvaluatorBase<Vec,Mat>(model),
    _queso_env(queso_env),
    _interface(argc,argv,model_comm,forward_run_input)
  {}

  template<class Vec,class Mat>
  void FullModelEvaluator<Vec,Mat>::compute_values( const std::vector<double>& param_values,
                                                    std::vector<double>& values ) const
  {
    // For this problem, we always have just two values
    queso_assert_equal_to( values.size(), 2 );

    // Update parameters will automatically resize
    std::vector<double> gamma_CN_params, gamma_N_params;
    this->_model.update_parameters(param_values,
                                   gamma_CN_params,
                                   gamma_N_params);

    this->_interface.set_gamma_CN_params( gamma_CN_params );
    this->_interface.set_gamma_N_params( gamma_N_params );

    try
      {
        this->_interface.solve();
      }
    catch(...)
      {
        if( this->_queso_env.subRank() == 0 )
          {
            std::cerr << "Caught exception in solver!" << std::endl;
            std::cerr << "Parameter Values:" << std::endl;

            for( std::vector<double>::const_iterator it = param_values.begin();
                 it != param_values.end(); ++it )
              std::cerr << "value = " << *it << std::endl;
          }

        values[0] = 0.0;
        values[1] = 0.0;

        this->_interface.reset_initial_guess();
      }

    values[0] = this->_interface.computed_mass_loss();
    values[1] = this->_interface.computed_average_n();

    this->_interface.reset_initial_guess();
  }

  // Instantiate GSL version of this class
  template class FullModelEvaluator<QUESO::GslVector,QUESO::GslMatrix>;

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
