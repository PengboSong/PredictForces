#include "read_pdb.h"
#include "addstr.h"

size_t read_resid(std::string line)
{
	return size_t(ssliceti(22, 26, line));
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
	res.resid = ssliceti(22, 26, line);
	res.x = sslicetd(30, 38, line);
	res.y = sslicetd(38, 46, line);
	res.z = sslicetd(46, 54, line);
	res.bfactor = sslicetd(60, 66, line);
	return res;
}

AtomInfo read_atom(std::string line)
{
	AtomInfo atom;
	atom.atomid = ssliceti(6, 11, line);
	atom.atomname = sslice(12, 16, line);
	atom.resname = sslice(17, 20, line);
	atom.chain = sslice(21, 22, line);
	atom.resid = ssliceti(22, 26, line);
	atom.x = sslicetd(30, 38, line);
	atom.y = sslicetd(38, 46, line);
	atom.z = sslicetd(46, 54, line);
	atom.bfactor = sslicetd(60, 66, line);

	return atom;
}