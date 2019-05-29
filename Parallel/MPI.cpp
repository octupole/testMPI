/*
 * MPI.cpp
 *
 *  Created on: May 24, 2019
 *      Author: marchi
 */

#include "MPI.h"
#include <mpi.h>
#include <cmath>
#include <boost/range/irange.hpp>


namespace Parallel {

vector<int> MPI::findSize3D(int n){
	int m=floor(pow(n,1/3));
	vector<int> result{m,m,m};

	if(m*m*m == n){
		return result;
	} else{
		while(true){
			for( size_t o{0};o<DIM;o++){
				result[o]++;
				if(result[XX]*result[YY]*result[ZZ] == n) return result;
				if(result[XX]*result[YY]*result[ZZ] > n) {
					result[o]--;
					return result;
				}

			}
		}
	}
};

vector<int> MPI::findSize3D(const int n,vector<int> vn){
	vector<int> result(vn);
	while(true){
		for( size_t o{0};o<DIM;o++){
			result[o]--;
			if(result[XX]*result[YY]*result[ZZ] <= n) return result;
		}
	}
};

void MPI::setDims(const int n, double cut, double a, double b, double c){
	vector<int> my{static_cast<int>(a/cut),static_cast<int>(b/cut),static_cast<int>(c/cut)};
	if(my[XX]*my[YY]*my[ZZ] == n){
		nc=my;
	}else if(my[XX]*my[YY]*my[ZZ] < n){
		int n_old=n;
		nc=my;
	}else if(my[XX]*my[YY]*my[ZZ] > n){
		nc=this->findSize3D(n,my);
	}
}
mpi::communicator MPI::ResizeWorld(mpi::communicator worldx0){
	if(world0.size() != nc[XX]*nc[YY]*nc[ZZ]){
		mpi::group local = worldx0.group();
		boost::integer_range<int> r = boost::irange(nc[XX]*nc[YY]*nc[ZZ]-1, worldx0.size()-1);
		mpi::group subgroup = local.exclude(r.begin(), r.end());
		mpi::communicator others{worldx0, subgroup};
		return others;
	} else{
		return world0;
	}

}

void MPI::PrintInfo(){
	if(world0.size() != world.size()){
		cout << " MPI size has been changed for Cartesian decomposition." << endl;
		cout << "\tOld size: "<< world0.size()<<endl;
		cout << "\tNew size: "<< world.size()<<": " << nc[XX]<< " X " << nc[YY]<< "  X " << nc[ZZ]<< " " <<endl;

	} else{
		cout << " MPI size for Cartesian decomposition." << endl;
		cout << "\tMPI size: "<< world.size()<<": " << nc[XX]<< " X " << nc[YY]<< "  X " << nc[ZZ]<< " " <<endl;
	}
}
void MPI::CartInit(){
	Dims[XX]={nc[XX], true};
	Dims[YY]={nc[YY], true};
	Dims[ZZ]={nc[ZZ], true};
	Cartx=new mpi::cartesian_communicator(world, mpi::cartesian_topology(Dims));

	for (int r = 0; r < Cartx->size(); ++r) {
		Cartx->barrier();
		if (r == Cartx->rank()) {
			std::vector<int> c = Cartx->coordinates(r);
			std::cout << "rk :" << r << " coords: "
					<< c[0] << ' ' << c[1] << ' ' << c[2] << '\n';
		}
	}
}
MPI::MPI(double cut, double a, double b, double c){

	if(world0.rank() != 0)
		std::cout.setstate(std::ios::failbit);

	this->setDims(world0.size(),cut,a,b,c);
	world=this->ResizeWorld(world0);
	myWorld=(world);
}
MPI::MPI(){
	if(world0.rank() != 0)
		std::cout.setstate(std::ios::failbit);

	nc=this->findSize3D(world0.size());
	world=this->ResizeWorld(world0);
	if((world) == MPI_COMM_NULL){
		   MPI_Finalize();
		   exit(0);
	}else{
		myWorld=(world);
	}
	int nrank;
	MPI_Comm_rank (myWorld, &nrank);
	cout << world.size() <<" " << nrank <<endl;
}

MPI::~MPI() {
	delete Cartx;
}

} /* namespace System */
