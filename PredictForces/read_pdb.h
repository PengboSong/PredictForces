#pragma once
#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

using boost::algorithm::trim;
using boost::lexical_cast;

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

std::string sslice(size_t begin, size_t size, std::string in);

size_t read_resid(std::string line);

std::string read_record(std::string line);

std::string read_resname(std::string line);

std::string read_chain(std::string line);

std::string read_atomname(std::string line);

ResInfo read_res(std::string line);

AtomInfo read_atom(std::string line);
