#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include "Pro.h"
#include "ProAnalysis.h"

std::list<double> cutoff_pocket = { 2.5, 3.0, 3.5, 4.0, 4.5, 5.0 };

std::vector<std::string> get_exclude_res()
{
	std::string exclbuf = "";
	std::cout << "Enter residues to be excluded: ";
	std::cin >> exclbuf;
	boost::trim(exclbuf);
	std::vector<std::string> excl = {};
	boost::split(excl, exclbuf, boost::is_any_of("\t, "));
	return excl;
}

int main()
{
	boost::filesystem::path datasetpath = "", proname = "";
	std::cout << "Enter dataset path: ";
	std::cin >> datasetpath;
	std::cout << "Enter protein family name: ";
	std::cin >> proname;

	boost::filesystem::path datadir = datasetpath / proname;
	if (!boost::filesystem::is_directory(datadir))
	{
		std::cout << "[Error] Can not find directory " << datadir << std::endl;
		exit(1);
	}

	std::string aponame = "", bindingname = "", allosteryname = "", complexname = "";
	std::cout << "Enter PDB name for apo state:";
	std::cin >> aponame;
	std::cout << "Enter PDB name for binding state:";
	std::cin >> bindingname;
	std::cout << "Enter PDB name for allostery state:";
	std::cin >> allosteryname;
	std::cout << "Enter PDB name for complex state:";
	std::cin >> complexname;

	Pro Apo;
	Pro Binding;
	Pro Allostery;
	Pro Complex;

	if (aponame.empty())
	{
		std::cout << "[Error] Apo state must be loaded." << std::endl;
		exit(1);
	}
	else
	{
		boost::filesystem::path apopath = datadir / aponame;
		std::vector<std::string> apoexcl = get_exclude_res();
		Apo = Pro::Pro(apopath.string(), false, apoexcl);
	}

	if (!bindingname.empty())
	{
		boost::filesystem::path bindingpath = datadir / bindingname;
		std::vector<std::string> bindingexcl = get_exclude_res();
		Binding = Pro::Pro(bindingpath.string(), true, bindingexcl);
	}
	if (!allosteryname.empty())
	{
		boost::filesystem::path allosterypath = datadir / allosteryname;
		std::vector<std::string> allosteryexcl = get_exclude_res();
		Allostery = Pro::Pro(allosterypath.string(), true, allosteryexcl);
	}
	if (!complexname.empty())
	{
		boost::filesystem::path complexpath = datadir / complexname;
		std::vector<std::string> complexexcl = get_exclude_res();
		Complex = Pro::Pro(complexpath.string(), true), complexexcl;
	}

	ProAnalysis Cycle(Apo, Pro::Pro(), Allostery, Complex);
	
	Cycle.interactive();
	
	return 0;
}