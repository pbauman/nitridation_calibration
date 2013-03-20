//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$ 
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef NITCAL_GAMMA_CN_LIKELIHOOD_H
#define NITCAL_GAMMA_CN_LIKELIHOOD_H

#ifdef HAVE_QUESO

#include "likelihood_base.h"
#include "simulation_interface.h"

namespace NitridationCalibration
{

  template<class Vec,class Mat>
  class GammaCNLikelihood : public LikelihoodBase<Vec,Mat>
  {
  public:

    GammaCNLikelihood( int argc,
		       char** argv,
		       MPI_Comm mpi_comm,
		       const std::string& input_filename,
		       const char* prefix, 
		       const uqVectorSetClass<Vec,Mat>& domain_set,
		       const bool returns_ln );

    virtual ~GammaCNLikelihood();

    virtual double evaluate_likelihood() const;

    virtual void update_parameters( const std::vector<double>& params ) const;

  protected:

    SimulationInterface _interface;

    double _sigma;
    double _sigma_sq;

    double _constant;

    double _data_mass_loss;
    
    double _gamma_nom;

  };

} // end namespace NitridationCalibration

#endif // HAVE_QUESO

#endif // NITCAL_GAMMA_CN_LIKELIHOOD_H
