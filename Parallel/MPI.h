/*
 * MPI.h
 *
 *  Created on: May 24, 2019
 *      Author: marchi
 */

#ifndef PARALLEL_MPI_H_
#define PARALLEL_MPI_H_
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
namespace mpi = boost::mpi;

namespace Parallel {

class MPI {
	mpi::communicator * worldx{nullptr};
	mpi::environment * envx{nullptr};
public:
	MPI();
	mpi::communicator & gWorld(){return *worldx;}
	virtual ~MPI();
};

} /* namespace System */

#endif /* PARALLEL_MPI_H_ */
