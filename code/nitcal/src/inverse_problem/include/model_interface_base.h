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

    ModelInterfaceBase( int argc, char** argv,
                        const QUESO::BaseEnvironment& queso_env,
                        const GetPot& forward_run_input );

    virtual ~ModelInterfaceBase(){};

    void solve( const std::vector<double>& param_values,
                std::vector<double>& values );

    const QUESO::VectorSpace<Vec,Mat>& param_space() const
    { return *(_param_space); }

    const QUESO::BoxSubset<Vec,Mat>& param_domain() const
    { return *(_param_domain.get()); }

  protected:

    virtual void update_parameters( const std::vector<double>& param_values ) =0;

    const QUESO::BaseEnvironment& _queso_env;

    SimulationInterface _interface;

    boost::scoped_ptr<QUESO::VectorSpace<Vec,Mat> > _param_space;
    boost::scoped_ptr<QUESO::BoxSubset<Vec,Mat> > _param_domain;

  private:

    ModelInterfaceBase();

  };

}// end namespace NitridationCalibration

#endif // NITCAL_MODEL_INTERFACE_BASE_H

#endif // NITCAL_HAVE_QUESO
