//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#ifndef NITCAL_MODEL_INTERPOLATION_BUILDER_H
#define NITCAL_MODEL_INTERPOLATION_BUILDER_H

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// NitCal
#include "model_evaluator_base.h"

// QUESO
#include "queso/InterpolationSurrogateBuilder.h"

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  class ModelInterpolationBuilder : public QUESO::InterpolationSurrogateBuilder<Vec,Mat>
  {
  public:

    ModelInterpolationBuilder( QUESO::InterpolationSurrogateDataSet<Vec,Mat>& data,
                               const ModelEvaluatorBase<Vec,Mat>& model );

    virtual ~ModelInterpolationBuilder(){};

    virtual void evaluate_model( const Vec& domainVector,
                                 std::vector<double>& values );

  private:

    const ModelEvaluatorBase<Vec,Mat>& _model;

    unsigned int _count;

  };

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO

#endif // NITCAL_MODEL_INTERPOLATION_BUILDER_H
