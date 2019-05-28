/*
 * MPI.cpp
 *
 *  Created on: May 24, 2019
 *      Author: marchi
 */

#include "MPI.h"
#include <mpi.h>

namespace Parallel {

MPI::MPI() {
	envx=new mpi::environment;
	worldx=new mpi::communicator;
	if(worldx->rank() != 0)
		std::cout.setstate(std::ios::failbit);
	// to return to active state: cout.clear()
	int world_size;
	MPI_Comm h=(*worldx);


}

MPI::~MPI() {
	delete envx;
	delete worldx;
}

} /* namespace System */
