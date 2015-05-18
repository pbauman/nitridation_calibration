//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

#ifndef NITCAL_SURROGATE_MODEL_COMPOSITION_H
#define NITCAL_SURROGATE_MODEL_COMPOSITION_H

// C++
#include <vector>

// Boost
#include "boost/scoped_ptr.hpp"

// NitCal
#include "model_interface_base.h"

// Foward declarations
class GetPot;

namespace QUESO
{
  class BaseEnvironment;

  template<class V, class M>
  class InterpolationSurrogateIOASCII;

  template<class V, class M>
  class LinearLagrangeInterpolationSurrogate;
}

namespace NitridationCalibration
{
  template<class Vec,class Mat>
  class SurrogateModelComposition
  {
  public:

    SurrogateModelComposition( const QUESO::BaseEnvironment& queso_env,
                               const GetPot& model_input );

    ~SurrogateModelComposition();

    void compute_values( const Vec& param_values,
                         Vec& model_output ) const;

    const ModelInterfaceBase<Vec,Mat>& get_model() const
    { return *(this->_model.get()); }

    const Vec& get_observations() const
    { return *(this->_observations.get()); }

    const Vec& get_covariance() const
    { return *(this->_covariance.get()); }

  protected:

    boost::scoped_ptr<ModelInterfaceBase<Vec,Mat> > _model;

    std::vector<GetPot*> _forward_model_inputs;

    boost::scoped_ptr<Vec> _observations;

    boost::scoped_ptr<Vec> _covariance;

    std::vector<QUESO::InterpolationSurrogateIOASCII<Vec,Mat>*> _interp_io;

    std::vector<QUESO::LinearLagrangeInterpolationSurrogate<Vec,Mat>*> _interp_surrogate;

  private:

    SurrogateModelComposition();

  };

} // end namespace NitridationCalibration

#endif // NITCAL_SURROGATE_MODEL_LIKELIHOOD_H

#endif // NITCAL_HAVE_QUESO
