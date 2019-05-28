//============================================================================
// Name        : testMPI.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <MPI.h>

using namespace std;

int main() {
//	Parallel::MPI my(10.0,100,100,100);
	Parallel::MPI my;
//	my.PrintInfo();
	my.CartInit();
	return 0;
}
