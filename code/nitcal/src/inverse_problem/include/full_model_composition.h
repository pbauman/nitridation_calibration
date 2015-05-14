//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

#ifndef NITCAL_FULL_MODEL_COMPOSITION_H
#define NITCAL_FULL_MODEL_COMPOSITION_H

#include "likelihood_comm_handler.h"
#include "model_interface_base.h"
#include "model_evaluator_base.h"

// Foward declarations
class GetPot;

namespace QUESO
{
  class BaseEnvironment;
}

namespace NitridationCalibration
{
  template<class Vec,class Mat>
  class FullModelComposition
  {
  public:

    FullModelComposition( int argc, char** argv,
                          const QUESO::BaseEnvironment& queso_env,
                          const GetPot& model_input );

    ~FullModelComposition();

  protected:

    boost::scoped_ptr<ModelInterfaceBase<Vec,Mat> > _model;

    LikelihoodCommHandler _comm_handler;

    boost::scoped_ptr<GetPot> _forward_run_input;

    boost::scoped_ptr<ModelEvaluatorBase<Vec,Mat> > _model_evaluator;

  private:

    FullModelComposition();

  };

} // end namespace NitridationCalibration

#endif // NITCAL_FULL_MODEL_LIKELIHOOD_H

#endif // NITCAL_HAVE_QUESO
