//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

// This class
#include "average_N_helper.h"

// GRINS
#include "grins/math_constants.h"

// libMesh
#include "libmesh/getpot.h"

namespace NitridationCalibration
{

  AverageNHelper::AverageNHelper( const std::string& input_filename )
  {
    GetPot input( input_filename );
    
    if( !input.have_variable( "AverageNLikelihood/sigma" ) )
      {
	std::cerr << "Error: Must specify sigma for mass loss likelihood."
		  << std::endl;
	libmesh_error();
      }

    const libMesh::Real sigma = input( "AverageNLikelihood/sigma", 0.0 );

    _sigma_sq = sigma*sigma;

    _constant = std::log(std::sqrt(GRINS::Constants::two_pi)*sigma);

    if( !input.have_variable( "AverageNLikelihood/data_value" ) )
      {
	std::cerr << "Error: Must specify data_value for mass loss likelihood."
		  << std::endl;
	libmesh_error();
      }

    _data = input( "AverageNLikelihood/data_value", 0.0 );

    return;
  }

  AverageNHelper::~AverageNHelper()
  {
    return;
  }

  libMesh::Real AverageNHelper::likelihood_value( const libMesh::Real computed_N_moles ) const
  {
    libMesh::Real likelihood_value = 0.0;

    const double tmp = computed_N_moles - _data;

    const double value = tmp*tmp/_sigma_sq;

    likelihood_value += value;

    likelihood_value += _constant;

    likelihood_value =  -likelihood_value/2.0;

    return likelihood_value;
  } 

} // end namespace NitridationCalibration
