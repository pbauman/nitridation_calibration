//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#ifndef NITCAL_FULL_MODEL_EVALUATOR_H
#define NITCAL_FULL_MODEL_EVALUATOR_H

#include "model_interface_base.h"
#include "model_evaluator_base.h"

// QUESO
#include <queso/VectorSet.h>
#include <queso/BoxSubset.h>

namespace NitridationCalibration
{
  template<class Vec,class Mat>
  class FullModelEvaluator : public ModelEvaluatorBase<Vec,Mat>
  {
  public:

    FullModelEvaluator( int argc, char** argv,
                        const QUESO::BaseEnvironment& queso_env,
                        const GetPot& forward_run_input,
                        MPI_Comm model_comm,
                        const ModelInterfaceBase<Vec,Mat>& model);

    virtual ~FullModelEvaluator(){};

    //! Compute output values from the model
    /*! This is const because we'll eventually put this in a likelihood
        which will have a const evaluate_model method */
    virtual void compute_values( const std::vector<double>& param_values,
                                 std::vector<double>& values ) const;

  protected:

    const QUESO::BaseEnvironment& _queso_env;

    SimulationInterface _interface;

  private:

    FullModelEvaluator();

  };

}// end namespace NitridationCalibration

#endif // NITCAL_FULL_MODEL_EVALUATOR_H
