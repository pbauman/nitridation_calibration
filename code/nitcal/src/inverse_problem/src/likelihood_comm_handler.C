//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

// This class
#include "likelihood_comm_handler.h"

namespace NitridationCalibration
{

  LikelihoodCommHandler::LikelihoodCommHandler( MPI_Comm& inter_chain_comm,
                                                const int n_datasets )
    : _inter_chain_comm(inter_chain_comm)
  {
    int inter_chain_size;
    MPI_Comm_size( _inter_chain_comm, &inter_chain_size );
    
    _n_procs_per_dataset = inter_chain_size/n_datasets;

    if( inter_chain_size % _n_procs_per_dataset != 0 )
      {
        std::cerr << "Error: The size of the Markov chain is not an integer multiple of"
                  << std::endl
                  << "       of the number of processes per forward problem."
                  << std::endl
                  << " inter_chain_size = " << inter_chain_size
                  << std::endl
                  << " n_procs_per_dataset = " << _n_procs_per_dataset
                  << std::endl;
        libmesh_error();
      }

    this->split_inter_chain_comm(n_datasets);

    return;
  }

  LikelihoodCommHandler::~LikelihoodCommHandler()
  {}

  MPI_Comm LikelihoodCommHandler::get_inter_chain_comm() const
  {
    return _inter_chain_comm;
  }

  MPI_Comm LikelihoodCommHandler::get_split_chain_comm() const
  {
    return _split_chain_comm;
  }

  MPI_Comm LikelihoodCommHandler::get_inter_chain_0_comm() const
  {
    return _inter_chain_0_comm;
  }

  int LikelihoodCommHandler::get_n_procs_per_dataset() const
  {
    return _n_procs_per_dataset;
  }

  void LikelihoodCommHandler::split_inter_chain_comm( const int n_datasets )
  {
    int inter_chain_rank;
    MPI_Comm_rank( _inter_chain_comm, &inter_chain_rank );
    
    // Processes with the same split_key will be in the same group
    // So, if we want more than 1 proc to run for the data sets
    // group them together.
    int split_key = inter_chain_rank/ _n_procs_per_dataset;

    int ierr = MPI_Comm_split( _inter_chain_comm, // Old communicator
                               split_key,
                               inter_chain_rank, // rank in old communicator
                               &_split_chain_comm ); // new communicator
    
    if( ierr != MPI_SUCCESS )
      {
        std::cerr << "Error occurred during MPI_Comm_split." << std::endl;
        libmesh_error();
      }

    // Now build up inter0 communicator for communicating between the masters of these
    // groups of processes
    std::vector<int> inter_chain_0_ranks;
    for( int i = 0; i < n_datasets;
         i++)
      {
        inter_chain_0_ranks.push_back( i*_n_procs_per_dataset);
      }

    MPI_Group  inter_chain_full_group;
    ierr = MPI_Comm_group( _inter_chain_comm, &inter_chain_full_group );

    ierr = MPI_Group_incl( inter_chain_full_group, 
                           n_datasets,
                           &inter_chain_0_ranks[0],
                           &_inter_chain_0_group );

    ierr = MPI_Comm_create( _inter_chain_comm,
                            _inter_chain_0_group,
                            &_inter_chain_0_comm );

    return;
  }


} // end namespace NitridationCalibration
