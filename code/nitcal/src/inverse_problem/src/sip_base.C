//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// This class
#include "sip_base.h"

// QUESO
#include "uqGslVector.h"
#include "uqGslMatrix.h"

namespace NitridationCalibration
{
  template<class Vec,class Mat>
  StatisticalInverseProblemBase<Vec,Mat>::StatisticalInverseProblemBase( uqBaseEnvironmentClass* env,
                                                                         const std::string& method,
                                                                         int argc,
                                                                         char** argv,
                                                                         const std::string& sip_input_filename )
    : QuesoStatisticalInverseProblemInterface<Vec,Mat>(env,method)
  {
    _sip_input.reset( new GetPot(sip_input_filename) );

    unsigned int n_datasets = (*_sip_input).vector_variable_size( "InverseProblem/datasets" );

    this->_comm_handler.reset( new LikelihoodCommHandler( env->subComm().Comm(), n_datasets ) );
    
    std::vector<std::string> datasets(n_datasets);
    for( unsigned int d = 0; d < n_datasets; d++ )
      {
        datasets[d] = (*_sip_input)( "InverseProblem/datasets", "DIE!", d );
      }

    int dataset_index = this->_comm_handler->get_dataset_index();

    _forward_run_input.reset( new GetPot(datasets[dataset_index]) );

    return;
  }


  // Instantiate GSL version of this class
  template class StatisticalInverseProblemBase<uqGslVectorClass,uqGslMatrixClass>;

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO
