//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#ifndef NITCAL_MODEL_EVALUATOR_BASE_H
#define NITCAL_MODEL_EVALUATOR_BASE_H

#include "model_interface_base.h"

// QUESO
#include <queso/VectorSet.h>
#include <queso/BoxSubset.h>

namespace NitridationCalibration
{
  template<class Vec,class Mat>
  class ModelEvaluatorBase
  {
  public:

    ModelEvaluatorBase( const ModelInterfaceBase<Vec,Mat>& model )
      : _model(model)
    {};

    virtual ~ModelEvaluatorBase(){};

    //! Compute output values from the model
    /*! This is const because we'll eventually put this in a likelihood
        which will have a const evaluate_model method */
    virtual void compute_values( const std::vector<double>& param_values,
                                 std::vector<double>& values ) const  =0;

  protected:

    const ModelInterfaceBase<Vec,Mat>& _model;

  private:

    ModelEvaluatorBase();

  };

}// end namespace NitridationCalibration

#endif // NITCAL_MODEL_EVALUATOR_BASE_H
