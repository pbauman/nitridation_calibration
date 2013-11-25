//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$ 
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef NITCAL_CONSTANT_GAMMA_CN_LIKELIHOOD_H
#define NITCAL_CONSTANT_GAMMA_CN_LIKELIHOOD_H

#include "nitcal_config.h"

#ifdef NITCAL_HAVE_QUESO

#include "likelihood_base.h"
#include "simulation_interface.h"
#include "mass_loss_helper.h"

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  class ConstantGammaCNLikelihood : public LikelihoodBase<Vec,Mat>
  {
  public:

    ConstantGammaCNLikelihood( int argc,
                               char** argv,
                               MPI_Comm mpi_comm,
                               const std::string& input_filename,
                               const char* prefix, 
                               const uqVectorSetClass<Vec,Mat>& domain_set,
                               const bool returns_ln );

    virtual ~ConstantGammaCNLikelihood();

    virtual double evaluate_likelihood() const;

    virtual void update_parameters( const std::vector<double>& params ) const;

  protected:

    MassLossHelper _mass_loss;
    
    double _gamma_nom;

  };

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO

#endif // NITCAL_CONSTANT_GAMMA_CN_LIKELIHOOD_H