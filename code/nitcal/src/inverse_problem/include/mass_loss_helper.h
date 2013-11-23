//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#ifndef NITCAL_MASS_LOSS_HELPER_H
#define NITCAL_MASS_LOSS_HELPER_H

// C++
#include <string>

// libMesh
#include "libmesh/libmesh_common.h"

namespace NitridationCalibration
{

  class MassLossHelper
  {
  public:

    MassLossHelper( const std::string& input_filename );
    ~MassLossHelper();

    libMesh::Real likelihood_value( const libMesh::Real computed_mass_loss ) const;

  private:

    MassLossHelper();

    libMesh::Real _sigma_sq;

    libMesh::Real _constant;

    libMesh::Real _data;

  };

} // end namespace NitridationCalibration

#endif // NITCAL_MASS_LOSS_HELPER_H
