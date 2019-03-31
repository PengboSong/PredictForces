#include <filesystem>

#include "Pro.h"
#include "ProAnalysis.h"

std::list<double> cutoff_pocket = { 2.5, 3.0, 3.5, 4.0, 4.5, 5.0 };

int main(int argc, char** argv)
{
	std::string tmp;
	std::string dataset = "C:\\Users\\china\\source\\repos\\PredictForces\\dataset\\";
	std::string metjapo = dataset + "MetJ\\1cmb.pdb";
	std::string metjc = dataset + "MetJ\\1cma.pdb";
	std::string metja = dataset + "MetJ\\1cmc.pdb";

	std::string p38apo = dataset + "p38\\1r39.pdb";
	std::string p38a1 = dataset + "p38\\1kv1.pdb";
	std::string p38a2 = dataset + "p38\\1kv2.pdb";

	Pro Apo(metjapo, false);
	Pro Alle(metja, true);
	Pro Complex(metjc, true);

	ProAnalysis Metj(Apo, Pro::Pro(), Alle, Complex);

	std::cout << "Using cutoff " << 4.0 << std::endl;
	Metj.set_learning_step(1e-4);
	Metj.set_convergence(0.1);
	Metj.gen_pocketA(2.0);
	Metj.test_pocketA();
	Metj.set_learning_step(1e-4);
	Metj.set_convergence(0.1);
	Metj.gen_pocketAS(2.6);
	Metj.test_pocketAS();

	/*
	for (std::list<double>::iterator it = cutoff_pocket.begin(); it != cutoff_pocket.end(); it++)
	{
		std::cout << "Using cutoff " << *it << std::endl;
		Metj.gen_pocketA(*it);
		Metj.test_pocketA();
	}

	for (std::list<double>::iterator it = cutoff_pocket.begin(); it != cutoff_pocket.end(); it++)
	{
		std::cout << "Using cutoff " << *it << std::endl;
		Metj.gen_pocketAS(*it);
		Metj.test_pocketAS();
	}*/

	// Metj.gen_pocketA(15.0);
	// Metj.gen_pocketAS(15.0);
	// Metj.gen_free_energy();
	
	std::cin >> tmp;
	return 0;
}