/********************************************************************
	created:	2012/11/10
	created:	10:11:2012   15:14
	filename: 	ParamTest.cpp
	file path:	GAG\test
	file base:	ParamTest
	file ext:	cpp
	author:		Han Hu
	
	purpose:	
*********************************************************************/

#include "GAGPL/MISC/Param.h"
#include <iostream>

using namespace param;


void printParam()
{
	// Print parameters.
	Param& param = Param::Instance();
	int max_charge = param.getParameter<int>("max_charge").first;
	int min_charge = param.getParameter<int>("min_charge").first;
	double tol = param.getParameter<double>("tolerance").first;
	int notexist = param.getParameter<int>("nothing").first;

	std::cout << "Max: " << max_charge << std::endl;
	std::cout << "Min: " << min_charge << std::endl;
	std::cout << "Tolerance: " << tol << std::endl;
	std::cout << "Nothing: " << notexist << std::endl;	
}

int main(int argc, char* argv[])
{
	

	// Set parameters.
	Param& param = Param::Instance();
	std::cout << "Get default value: " << std::endl;
	printParam();

	std::cout << "Set parameters: " << std::endl;
	param.setParameter("max_charge", 3);
	param.setParameter("min_charge", 1);
	param.setParameter<double>("tolerance", 5e-6);
	param.setParameter<double>("A2tolerance", 1e-5);

	printParam();

	// TBD
	// Read parameters from file.

	// Change parameters.
	std::cout << "Change parameters:" << std::endl;
	param.setParameter("max_charge", 4);
	//max_charge = param.getParameter<int>("max_charge").first;
	//std::cout << "New Max: " << max_charge << std::endl;
	printParam();
}