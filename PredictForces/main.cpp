#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include "Pro.h"
#include "ProAnalysis.h"

using boost::algorithm::trim;

std::list<double> cutoff_pocket = { 2.5, 3.0, 3.5, 4.0, 4.5, 5.0 };

std::vector<std::string> get_exclude_res()
{
	std::string exclbuf = "";
	std::cout << "Enter residues to be excluded: ";
	std::getline(std::cin, exclbuf);
	trim(exclbuf);
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

	Pro Apo;
	Pro Binding;
	Pro Allostery;
	Pro Complex;

	std::cin.get();
	std::string aponame = "", bindingname = "", allosteryname = "", complexname = "";
	std::cout << "Enter PDB name for apo state:";
	std::getline(std::cin, aponame);
	trim(aponame);

	if (aponame.empty())
	{
		std::cout << "[Error] Apo state must be loaded." << std::endl;
		exit(1);
	}
	else
	{
		if (!boost::algorithm::ends_with(aponame, ".pdb"))
			aponame += ".pdb";
		boost::filesystem::path apopath = datadir / aponame;
		std::vector<std::string> apoexcl = get_exclude_res();
		Apo = Pro(apopath.string(), false, apoexcl);
	}

	std::cout << "Enter PDB name for binding state:";
	std::getline(std::cin, bindingname);
	trim(bindingname);
	if (!bindingname.empty())
	{
		if (!boost::algorithm::ends_with(bindingname, ".pdb"))
			bindingname += ".pdb";
		boost::filesystem::path bindingpath = datadir / bindingname;
		std::vector<std::string> bindingexcl = get_exclude_res();
		Binding = Pro(bindingpath.string(), true, bindingexcl);
	}

	std::cout << "Enter PDB name for allostery state:";
	std::getline(std::cin, allosteryname);
	trim(allosteryname);
	if (!allosteryname.empty())
	{
		if (!boost::algorithm::ends_with(allosteryname, ".pdb"))
			allosteryname += ".pdb";
		boost::filesystem::path allosterypath = datadir / allosteryname;
		std::vector<std::string> allosteryexcl = get_exclude_res();
		Allostery = Pro(allosterypath.string(), true, allosteryexcl);
	}

	std::cout << "Enter PDB name for complex state:";
	std::getline(std::cin, complexname);
	trim(complexname);
	if (!complexname.empty())
	{
		if (!boost::algorithm::ends_with(complexname, ".pdb"))
			complexname += ".pdb";
		boost::filesystem::path complexpath = datadir / complexname;
		std::vector<std::string> complexexcl = get_exclude_res();
		Complex = Pro(complexpath.string(), true, complexexcl);
	}

	ProAnalysis Cycle(Apo, Binding, Allostery, Complex);
	
	Cycle.interactive();
	
	return 0;
}