//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

// This class
#include "model_interface_base.h"

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// QUESO
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  ModelInterfaceBase<Vec,Mat>::ModelInterfaceBase( int argc, char** argv,
                                                   const QUESO::BaseEnvironment& queso_env,
                                                   const GetPot& forward_run_input )
    : _queso_env(queso_env),
      _interface(argc,argv,queso_env.subComm().Comm(),forward_run_input),
      _param_space(NULL),
      _param_domain(NULL)
  {}

  template<class Vec,class Mat>
  void ModelInterfaceBase<Vec,Mat>::solve( const std::vector<double>& param_values,
                                           std::vector<double>& values )
  {
    queso_assert_equal_to( values.size(), 2 );

    this->update_parameters(param_values);

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
  template class ModelInterfaceBase<QUESO::GslVector,QUESO::GslMatrix>;

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
