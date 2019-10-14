#include "read_pdb.h"

std::string sslice(size_t begin, size_t end, std::string in)
{
	std::string out = in.substr(begin, end - begin);
	boost::algorithm::trim(out);
	return out;
}

size_t read_resid(std::string line)
{
	return boost::lexical_cast<size_t>(sslice(22, 26, line));
}

std::string read_record(std::string line)
{
	return sslice(0, 6, line);
}

std::string read_resname(std::string line)
{
	return sslice(17, 20, line);
}

std::string read_chain(std::string line)
{
	return sslice(21, 22, line);
}

std::string read_atomname(std::string line)
{
	return sslice(12, 16, line);
}

ResInfo read_res(std::string line)
{
	ResInfo res;
	res.resname = sslice(17, 20, line);
	res.chain = sslice(21, 22, line);
	res.resid = boost::lexical_cast<size_t>(sslice(22, 26, line));
	res.x = boost::lexical_cast<double>(sslice(30, 38, line));
	res.y = boost::lexical_cast<double>(sslice(38, 46, line));
	res.z = boost::lexical_cast<double>(sslice(46, 54, line));
	res.bfactor = boost::lexical_cast<double>(sslice(60, 66, line));

	return res;
}

AtomInfo read_atom(std::string line)
{
	AtomInfo atom;
	atom.resid = boost::lexical_cast<size_t>(sslice(6, 11, line));
	atom.atomname = sslice(12, 16, line);
	atom.resname = sslice(17, 20, line);
	atom.chain = sslice(21, 22, line);
	atom.resid = boost::lexical_cast<size_t>(sslice(22, 26, line));
	atom.x = boost::lexical_cast<double>(sslice(30, 38, line));
	atom.y = boost::lexical_cast<double>(sslice(38, 46, line));
	atom.z = boost::lexical_cast<double>(sslice(46, 54, line));
	atom.bfactor = boost::lexical_cast<double>(sslice(60, 66, line));

	return atom;
}