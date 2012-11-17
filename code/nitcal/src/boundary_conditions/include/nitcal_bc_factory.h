//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef NITCAL_BC_FACTORY_H
#define NITCAL_BC_FACTORY_H

#include <string>

// libMesh
#include "getpot.h"
#include "equation_systems.h"
#include "system.h"
#include "zero_function.h"

// GRINS
#include "bc_types.h"
#include "bc_factory.h"

// NitridationCalibration
#include "nitcal_catalytic_wall.h"
#include "tube_twall.h"

namespace NitridationCalibration
{
  //! Object for constructing NitridationCalibration specific boundary condition function objects.
  class BoundaryConditionsFactory : public GRINS::BoundaryConditionsFactory
  {
  public:
    
    BoundaryConditionsFactory( const GetPot& input );

    virtual ~BoundaryConditionsFactory();
    
    //! Builds the NitridationCalibration::
    virtual std::multimap< GRINS::PhysicsName, GRINS::DBCContainer > build_dirichlet( );

    //! Builds the NitridationCalibration::
    virtual std::map< GRINS::PhysicsName, GRINS::NBCContainer > build_neumann( libMesh::EquationSystems& es );

  }; // class BoundaryConditionsFactory

} // namespace NitridationCalibration

#endif //NITCAL_BC_FACTORY_H
