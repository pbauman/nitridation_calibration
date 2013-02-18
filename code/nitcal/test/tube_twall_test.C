//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// NitCal - Nitridation Calibration 
//
// Copyright (C) 2012 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include <iomanip>
#include "tube_twall.h"

int main()
{
  GetPot input( "./input_files/tube_twall_test.in" );
  
  NitridationCalibration::TubeTempBC tube(input);

  libMesh::Point p0( 0.10 );
  libMesh::Point p1( 0.30 );
  libMesh::Point p2( 0.70 );
  libMesh::Point p3( 0.80 );
  
  const double tol = 1.0e-15;

  int return_flag = 0;

  {

    double T0 = tube( p0 );
    double T1 = tube( p1 );
    double T2 = tube( p2 );
    double T3 = tube( p3 );
  
    /*
      std::cout << std::setprecision(16) << std::scientific
      << "T0 = " << T0 << std::endl
      << "T1 = " << T1 << std::endl
      << "T2 = " << T2 << std::endl
      << "T3 = " << T3 << std::endl;
    */

    const double T0_reg = 3.1931147540983608e+02;
    const double T1_reg = 9.9609999999999980e+02;
    const double T2_reg = 6.2117142857142869e+02;
    const double T3_reg = 4.5390624999999994e+02;

    const double tol = 1.0e-15;

    if( std::fabs( (T0_reg - T0)/T0_reg ) > tol )
      {
	std::cerr << "Error: T0 mismatch." << std::endl
		  << "T0 = " << T0 << std::endl
		  << "T0_reg = " << T0_reg << std::endl;
	return_flag = 1;
      }

    if( std::fabs( (T1_reg - T1)/T1_reg ) > tol )
      {
	std::cerr << "Error: T1 mismatch." << std::endl
		  << "T1 = " << T1 << std::endl
		  << "T1_reg = " << T1_reg << std::endl;
	return_flag = 1;
      }

    if( std::fabs( (T2_reg - T2)/T2_reg ) > tol )
      {
	std::cerr << "Error: T2 mismatch." << std::endl
		  << "T2 = " << T2 << std::endl
		  << "T2_reg = " << T2_reg << std::endl;
	return_flag = 1;
      }

    if( std::fabs( (T3_reg - T3)/T3_reg ) > tol )
      {
	std::cerr << "Error: T3 mismatch." << std::endl
		  << "T3 = " << T3 << std::endl
		  << "T3_reg = " << T3_reg << std::endl;
	return_flag = 1;
      }

  }

  libMesh::AutoPtr<libMesh::FunctionBase<libMesh::Real> > tube_clone = tube.clone();

  {

    double T0 = (*tube_clone)( p0 );
    double T1 = (*tube_clone)( p1 );
    double T2 = (*tube_clone)( p2 );
    double T3 = (*tube_clone)( p3 );
  
    /*
      std::cout << std::setprecision(16) << std::scientific
      << "T0 = " << T0 << std::endl
      << "T1 = " << T1 << std::endl
      << "T2 = " << T2 << std::endl
      << "T3 = " << T3 << std::endl;
    */

    const double T0_reg = 3.1931147540983608e+02;
    const double T1_reg = 9.9609999999999980e+02;
    const double T2_reg = 6.2117142857142869e+02;
    const double T3_reg = 4.5390624999999994e+02;

    if( std::fabs( (T0_reg - T0)/T0_reg ) > tol )
      {
	std::cerr << "Error: T0 mismatch." << std::endl
		  << "T0 = " << T0 << std::endl
		  << "T0_reg = " << T0_reg << std::endl;
	return_flag = 1;
      }

    if( std::fabs( (T1_reg - T1)/T1_reg ) > tol )
      {
	std::cerr << "Error: T1 mismatch." << std::endl
		  << "T1 = " << T1 << std::endl
		  << "T1_reg = " << T1_reg << std::endl;
	return_flag = 1;
      }

    if( std::fabs( (T2_reg - T2)/T2_reg ) > tol )
      {
	std::cerr << "Error: T2 mismatch." << std::endl
		  << "T2 = " << T2 << std::endl
		  << "T2_reg = " << T2_reg << std::endl;
	return_flag = 1;
      }

    if( std::fabs( (T3_reg - T3)/T3_reg ) > tol )
      {
	std::cerr << "Error: T3 mismatch." << std::endl
		  << "T3 = " << T3 << std::endl
		  << "T3_reg = " << T3_reg << std::endl;
	return_flag = 1;
      }

  }

  return return_flag;
}
