//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// This class
#include "full_model_composition.h"

//NitCal
#include "full_model_evaluator.h"
#include "model_builder.h"

namespace NitridationCalibration
{
  template<class Vec,class Mat>
  FullModelComposition<Vec,Mat>::FullModelComposition( int argc, char** argv,
                                                       const QUESO::BaseEnvironment& queso_env,
                                                       const GetPot& model_input )
    : _model(ModelBuilder<Vec,Mat>::build_model(queso_env,model_input)),
      _comm_handler(queso_env.subComm().Comm(),
                    model_input.vector_variable_size("Likelihood/datasets") )
  {
    unsigned int n_datasets = model_input.vector_variable_size("Likelihood/datasets");

    std::vector<std::string> datasets(n_datasets);
    for( unsigned int d = 0; d < n_datasets; d++ )
      {
        datasets[d] = model_input( "Likelihood/datasets", "DIE!", d );
      }

    int dataset_index = this->_comm_handler->get_dataset_index();

    _forward_run_input.reset( new GetPot(datasets[dataset_index]) );

    _model_evaluator.reset( new FullModelEvaluator<Vec,Mat>(argc,argv,
                                                            queso_env,
                                                            *(_forward_run_input.get()),
                                                            _comm_handler.get_split_chain_comm(),
                                                            *(_model.get())) );

  }

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
