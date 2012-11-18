//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef NITCAL_CATALYTIC_WALL_H
#define NITCAL_CATALYTIC_WALL_H

namespace NitridationCalibration
{

  class CatalyticWall : public GRINS::NeumannFuncObj
  {
  public:

    CatalyticWall( libMesh::EquationSystems& es );
    virtual ~CatalyticWall();

    //! Returns the value of the implemented Neumann boundary condition
    /*! This will leverage the FEMContext to get variable values and derivatives through the
      side_value, side_gradient, etc. interfaces, for each quadrature point qp. */
    virtual libMesh::Point value( const libMesh::FEMContext& context, const unsigned int qp );

    //! Returns the value of the implemented Neumann boundary condition
    /*! This will leverage the FEMContext to get variable values and derivatives through the
      side_value, side_gradient, etc. interfaces, for each quadrature point qp. 
      Returns the normal component of the Neumann value. Only to be used when flux vector is
      formulated implicitly in terms of normal component. By default, does nothing since
      it's only applicable in special cases. */
    virtual Real normal_value( const libMesh::FEMContext& context, const unsigned int qp );

    //! Returns the derivative with respect to the primary variable of the implemented Neumann boundary condition.
    /*! This will leverage the FEMContext to get variable values and derivatives through the
      side_value, side_gradient, etc. interfaces, for each quadrature point qp. */
    virtual libMesh::Point derivative( const libMesh::FEMContext& context, const unsigned qp );

    //! Returns the derivative with respect to the primary variable of the implemented Neumann boundary condition.
    /*! This will leverage the FEMContext to get variable values and derivatives through the
      side_value, side_gradient, etc. interfaces, for each quadrature point qp.
      Returns the normal component of the Neumann value. Only to be used when flux vector is
      formulated implicitly in terms of normal component. By default, does nothing since
      it's only applicable in special cases. */
    virtual Real normal_derivative( const libMesh::FEMContext& context, const unsigned qp );


  protected:

    const GRINS::VariableIndex _T_var;
    const GRINS::VariableIndex _N_species_var;
    const GRINS::VariableIndex _CN_species_var;

    const GRINS::ChemicalMixture& _chem_mixture;
  };

}// namespace NitridationCalibration

#endif //NITCAL_CATALYTIC_WALL_H
