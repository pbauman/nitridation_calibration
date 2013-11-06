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

#ifndef NITCAL_QOI_FACTORY_H
#define NITCAL_QOI_FACTORY_H

// GRINS
#include "grins/qoi_factory.h"

namespace NitridationCalibration
{
  class QoIFactory : public GRINS::QoIFactory
  {
  public:

    QoIFactory();
    
    virtual ~QoIFactory();

  protected:

    virtual void add_qoi( const GetPot& input,
			  const std::string& qoi_name,
			  std::tr1::shared_ptr<GRINS::CompositeQoI>& qois );

  };

} // end namespace NitridationCalibration

#endif // NITCAL_QOI_FACTORY_H
