//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "nitcal_bc_factory.h"

namespace NitridationCalibration
{
  BoundaryConditionsFactory::BoundaryConditionsFactory( const GetPot& input )
    : _input(input)
  {
    return;
  }

  BoundaryConditionsFactory::~BoundaryConditionsFactory()
  {
    return;
  }

  std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > BoundaryConditionsFactory::build_dirichlet( )
  {
    std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > dbc;

    return dbc;
  }

  std::map< GRINS::PhysicsName, GRINS::NBCContainer > BoundaryConditionsFactory::build_neumann( libMesh::EquationSystems& es )
  {
    std::map< GRINS::PhysicsName, GRINS::NBCContainer > nbc;

    return nbc;
  }

} // NitridationCalibration
