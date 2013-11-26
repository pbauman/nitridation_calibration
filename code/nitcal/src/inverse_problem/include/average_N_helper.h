//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#ifndef NITCAL_AVERAGE_N_HELPER_H
#define NITCAL_AVERAGE_N_HELPER_H

// C++
#include <string>

// libMesh
#include "libmesh/libmesh_common.h"

class GetPot;

namespace NitridationCalibration
{

  class AverageNHelper
  {
  public:

    AverageNHelper( const GetPot& input );
    ~AverageNHelper();

    libMesh::Real likelihood_value( const libMesh::Real computed_mass_loss ) const;

  private:

    AverageNHelper();

    libMesh::Real _sigma_sq;

    libMesh::Real _constant;

    libMesh::Real _data;

  };

} // end namespace NitridationCalibration

#endif // NITCAL_AVERAGE_N_HELPER_H
