//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// NitCal - Nitridation Calibration 
//
// Copyright (C) 2012-2013 The PECOS Development Team
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

// This class
#include "qoi_factory.h"

// NitCal
#include "mass_loss.h"
#include "average_N_mole_fraction.h"

// GRINS
#include "grins/qoi_base.h"

namespace NitridationCalibration
{
  QoIFactory::QoIFactory()
    : GRINS::QoIFactory()
  {
    return;
  }

  QoIFactory::~QoIFactory()
  {
    return;
  }

  void QoIFactory::add_qoi( const GetPot& input,
			    const std::string& qoi_name,
			    std::tr1::shared_ptr<GRINS::CompositeQoI>& qois )
  {
    if( qoi_name == std::string("MassLoss") )
      {
        GRINS::QoIBase* qoi = new MassLoss( std::string("MassLoss") );
        qois->add_qoi(*qoi);
      }
    
    else if( qoi_name == std::string("AverageNMoleFraction") )
      {
        GRINS::QoIBase* qoi = new AverageNMoleFraction( std::string("AverageNMoleFraction") );
        qois->add_qoi(*qoi);
      }
    
    else
      {
	GRINS::QoIFactory::add_qoi( input, qoi_name, qois );
      }

    return;
  }

} // end namespace NitridationCalibration
