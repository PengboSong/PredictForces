#include "addstr.h"
#include "ReadPro.h"

ReadPro::ReadPro(std::string line)
{
	record = sslice(0, 6, line);
	if (record == "ATOM" || record == "HETATM")
	{
		atomid = ssliceti(6, 11, line);
		atomname = sslice(12, 16, line);
		resname = sslice(17, 20, line);
		chain = sslice(21, 22, line);
		resid = ssliceti(22, 26, line);
		x = sslicetd(30, 38, line);
		y = sslicetd(38, 46, line);
		z = sslicetd(46, 54, line);
		occup = sslicetd(54, 60, line);
		bfactor = sslicetd(60, 66, line);
		if (full_mode_flag)
		{
			location = sslice(16, 17, line);
			code = sslice(26, 27, line);
			segment = sslice(72, 76, line);
			element = sslice(76, 78, line);
		}
	}
}

ReadPro::~ReadPro()
{
}


void ReadPro::activate_full_mode()
{
	full_mode_flag = true;
}

void ReadPro::deactivate_full_mode()
{
	full_mode_flag = false;
}
