#pragma once
#include <string>

struct ResInfo
{
	std::string resname;
	std::string chain;
	size_t resid;
	double x;
	double y;
	double z;
	double bfactor;
};
struct AtomInfo
{
	std::string atomname;
	std::string resname;
	std::string chain;
	size_t atomid;
	size_t resid;
	double x;
	double y;
	double z;
	double bfactor;
};

size_t read_resid(std::string line);

std::string read_record(std::string line);

std::string read_resname(std::string line);

std::string read_chain(std::string line);

std::string read_atomname(std::string line);

ResInfo read_res(std::string line);

AtomInfo read_atom(std::string line);
