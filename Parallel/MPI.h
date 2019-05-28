/*
 * MPI.h
 *
 *  Created on: May 24, 2019
 *      Author: marchi
 */

#ifndef PARALLEL_MPI_H_
#define PARALLEL_MPI_H_
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/cartesian_communicator.hpp>
#include <vector>
#include <iostream>

using std::vector;
using std::cout;using std::endl;
namespace mpi = boost::mpi;

namespace Parallel {
const int DIM{3};
const int XX{0},YY{1},ZZ{2};
class MPI {
	mpi::communicator world0;
	mpi::communicator world;
	mpi::cartesian_communicator * Cartx;
	mpi::cartesian_dimension * Dims[3];
	mpi::environment env;
	MPI_Comm myWorld;
	vector<int> nc{0,0,0};
	void setDims(int, double, double, double, double);
	vector<int> findSize3D(const int, vector<int>);
	vector<int> findSize3D(const int);
	mpi::communicator ResizeWorld(mpi::communicator);
public:
	MPI();
	MPI(double,double,double,double);
	void CartInit();
	mpi::communicator & gWorld(){return world;}
	void PrintInfo();
	virtual ~MPI();
};

} /* namespace System */

#endif /* PARALLEL_MPI_H_ */
