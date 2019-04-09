#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include "manageIO.h"
#include "Pro.h"
#include "ProAnalysis.h"

using namespace std;
using namespace boost::algorithm;
using boost::filesystem::path;

list<double> cutoff_pocket = { 2.5, 3.0, 3.5, 4.0, 4.5, 5.0 };

vector<string> get_exclude_res()
{
	string exclbuf = "";
	cout << "Enter residues to be excluded: ";
	getline(cin, exclbuf);
	trim(exclbuf);
	vector<string> excl = {};
	boost::split(excl, exclbuf, boost::is_any_of("\t, "));
	return excl;
}

int main()
{
	path datasetpath = "", proname = "";
	cout << "Enter dataset path: ";
	cin >> datasetpath;
	cout << "Enter protein family name: ";
	cin >> proname;

	path datadir = datasetpath / proname;
	if (!boost::filesystem::is_directory(datadir))
		handle_error("Can not find directory " + datadir.string());

	Pro Apo;
	Pro Binding;
	Pro Allostery;
	Pro Complex;

	cin.get();
	string aponame = "", bindingname = "", allosteryname = "", complexname = "";
	cout << "Enter PDB name for apo state:";
	getline(cin, aponame);
	trim(aponame);

	if (aponame.empty())
		handle_error("Apo state must be loaded.");
	else
	{
		if (!ends_with(aponame, ".pdb"))
			aponame += ".pdb";
		path apopath = datadir / aponame;
		vector<string> apoexcl = get_exclude_res();
		Apo = Pro(apopath.string(), false, apoexcl);
	}

	cout << "Enter PDB name for binding state:";
	getline(cin, bindingname);
	trim(bindingname);
	if (!bindingname.empty())
	{
		if (!boost::algorithm::ends_with(bindingname, ".pdb"))
			bindingname += ".pdb";
		path bindingpath = datadir / bindingname;
		vector<string> bindingexcl = get_exclude_res();
		Binding = Pro(bindingpath.string(), true, bindingexcl);
	}

	cout << "Enter PDB name for allostery state:";
	getline(cin, allosteryname);
	trim(allosteryname);
	if (!allosteryname.empty())
	{
		if (!boost::algorithm::ends_with(allosteryname, ".pdb"))
			allosteryname += ".pdb";
		path allosterypath = datadir / allosteryname;
		vector<string> allosteryexcl = get_exclude_res();
		Allostery = Pro(allosterypath.string(), true, allosteryexcl);
	}

	cout << "Enter PDB name for complex state:";
	getline(cin, complexname);
	trim(complexname);
	if (!complexname.empty())
	{
		if (!boost::algorithm::ends_with(complexname, ".pdb"))
			complexname += ".pdb";
		path complexpath = datadir / complexname;
		vector<string> complexexcl = get_exclude_res();
		Complex = Pro(complexpath.string(), true, complexexcl);
	}

	ProAnalysis Cycle(Apo, Binding, Allostery, Complex);
	
	Cycle.interactive();
	
	return 0;
}