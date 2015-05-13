//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#ifndef NITCAL_MODEL_INTERFACE_BASE_H
#define NITCAL_MODEL_INTERFACE_BASE_H


#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

// Boost
#include "boost/scoped_ptr.hpp"

#include "simulation_interface.h"

// QUESO
#include <queso/VectorSet.h>
#include <queso/BoxSubset.h>

namespace NitridationCalibration
{
  template<class Vec,class Mat>
  class ModelInterfaceBase
  {
  public:

    ModelInterfaceBase()
      : _param_space(NULL),
        _param_domain(NULL)
    {}

    virtual ~ModelInterfaceBase(){};

    virtual void update_parameters( const std::vector<double>& param_values,
                                    std::vector<double>& gamma_CN_params,
                                    std::vector<double>& gamma_N_params ) const =0;

    const QUESO::VectorSpace<Vec,Mat>& param_space() const
    { return *(_param_space); }

    const QUESO::BoxSubset<Vec,Mat>& param_domain() const
    { return *(_param_domain.get()); }

  protected:

    boost::scoped_ptr<QUESO::VectorSpace<Vec,Mat> > _param_space;
    boost::scoped_ptr<QUESO::BoxSubset<Vec,Mat> > _param_domain;

  };

}// end namespace NitridationCalibration

#endif // NITCAL_MODEL_INTERFACE_BASE_H

#endif // NITCAL_HAVE_QUESO
