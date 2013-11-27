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
                               const GetPot& sip_input,
                               const GetPot& forward_run_input,
                               const LikelihoodCommHandler& comm_handler,
                               const char* prefix, 
                               const uqVectorSetClass<Vec,Mat>& domain_set,
                               const bool returns_ln );

    virtual ~ConstantGammaCNLikelihood();

    virtual double evaluate_likelihood() const;

    virtual void update_parameters( const std::vector<double>& params ) const;

  protected:

    MassLossHelper _mass_loss;
    
    double _gamma_nom;

  private:

    ConstantGammaCNLikelihood();

  };

} // end namespace NitridationCalibration

#endif // NITCAL_HAVE_QUESO

#endif // NITCAL_CONSTANT_GAMMA_CN_LIKELIHOOD_H
