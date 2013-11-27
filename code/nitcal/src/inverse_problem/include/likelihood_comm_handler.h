//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-

#ifndef NITCAL_LIKELHIHOOD_COMM_HANDLER_H
#define NITCAL_LIKELHIHOOD_COMM_HANDLER_H

#include "libmesh/libmesh.h"

namespace NitridationCalibration
{

  class LikelihoodCommHandler
  {
  public:

    LikelihoodCommHandler( MPI_Comm inter_chain_comm,
                           const int n_datasets  );
    ~LikelihoodCommHandler();

    MPI_Comm get_split_chain_comm() const;
    MPI_Comm get_inter_chain_comm() const;
    MPI_Comm get_inter_chain_0_comm() const;

    int get_n_procs_per_dataset() const;

    int get_dataset_index() const;

  protected:

    void split_inter_chain_comm( const int n_datasets );

    // This is the communicator for all procs within a markov chain
    // or a "QUESO subenvironment"
    MPI_Comm _inter_chain_comm;

    // This is the communicator for the split QUESO subenvironment
    // for running many experiments in parallel within likihood
    MPI_Comm _split_chain_comm;

    int _n_procs_per_dataset;

    MPI_Group _inter_chain_0_group;
    MPI_Comm  _inter_chain_0_comm;

    int _dataset_index;

  private:

    LikelihoodCommHandler();

  };

} // end namespace NitridationCalibration

#endif // NITCAL_LIKELHIHOOD_COMM_HANDLER_H
