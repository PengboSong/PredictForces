#pragma once
#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include "defines.h"

std::string sslice(size_t begin, size_t size, std::string in);

size_t read_resid(std::string line);

std::string read_record(std::string line);

std::string read_resname(std::string line);

std::string read_chain(std::string line);

std::string read_atomname(std::string line);

ResInfo read_res(std::string line);

AtomInfo read_atom(std::string line);
