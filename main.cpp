//============================================================================
// Name        : testMPI.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <vector>
#include "../Parallel/MPI.h"
#include <MPIconfig.hpp>
#include <map>

using namespace std;

int main() {

//	Parallel::MPI my(10.0,100,100,100);
	Parallel::MPI my;
//	my.PrintInfo();
	my.CartInit();
	map<int,int> outbuf;
	for(int o{0};o<Parallel::CARTDIRS;o++){
		auto Cart=my.gWorld();
		if(o%2 == 0) outbuf[o]=Cart.rank();
		else outbuf[o]=-Cart.rank();
	}

//	my.CartSend(outbuf);

//	my.PrintInfo();
	return 0;
}
