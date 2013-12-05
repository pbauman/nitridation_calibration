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

// Antioch
#include "antioch/vector_utils_decl.h"

// NitCal
#include "mass_loss_catalytic.h"
#include "average_N_mole_fraction.h"
#include "qoi_names.h"

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
    if( qoi_name == average_N_mole_fraction )
      {
        GRINS::QoIBase* qoi = new AverageNMoleFraction( average_N_mole_fraction );
        qois->add_qoi(*qoi);
      }
    
    else if( qoi_name == mass_loss_catalytic  )
      {
        GRINS::QoIBase* qoi = new MassLossCatalytic( mass_loss_catalytic );
        qois->add_qoi(*qoi);
      }

    else
      {
	GRINS::QoIFactory::add_qoi( input, qoi_name, qois );
      }

    return;
  }

} // end namespace NitridationCalibration
