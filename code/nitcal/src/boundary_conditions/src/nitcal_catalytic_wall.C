//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "nitcal_catalytic_wall.h"

namespace NitridationCalibration
{

  CatalyticWall::CatalyticWall( libMesh::EquationSystems& es )
  {
    const libMesh::System& system = es.get_system(_system_name);

    _T_var = system.variable_number( _T_var_name );
  

    _T_var = 

    return;
  }

  CatalyticWall::~CatalyticWall()
  {
    return;
  }

} //namespace NitridationCalibration
