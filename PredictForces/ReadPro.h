#pragma once

class ReadPro
{
public:
	ReadPro(std::string line);
	~ReadPro();

	void activate_full_mode();
	void deactivate_full_mode();

private:
	std::string record;
	size_t atomid;
	std::string atomname;
	std::string location;
	std::string resname;
	std::string chain;
	size_t resid;
	std::string code;
	double x;
	double y;
	double z;
	double occup;
	double bfactor;
	std::string segment;
	std::string element;

	bool is_pro = false;
	bool is_ligand = false;
	bool is_water = false;

	bool full_mode_flag = false;
};

