#pragma once
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <utility>
#include <limits>
#include <boost/format.hpp>
#include <Eigen/Dense>

#include "defines.h"
#include "handle_io.h"
#include "read_pdb.h"

using namespace Eigen;

class Pro
{
public:
	Pro();
	Pro(std::string fpath, ProConfigs configs);
	~Pro();

	bool empty() const { return resn == 0; }

	bool has_ligand() const { return with_ligand_flag; }

	bool has_res(size_t id) const { return id < resn; }

	std::string get_resname(size_t id) const
	{
		if (id < resn)
			return pro.at(id).resname;
		else
			return std::string();
	}

	std::string get_chain(size_t id) const
	{
		if (id < resn)
			return pro.at(id).chain;
		else
			return std::string();
	}

	size_t get_resid(size_t id) const
	{
		if (id < resn)
			return pro.at(id).resid;
		else
			return 0;
	}

	double get_x(size_t id) const
	{
		if (id < resn)
			return pro.at(id).x;
		else
			return 0.0;
	}

	double get_y(size_t id) const
	{
		if (id < resn)
			return pro.at(id).y;
		else
			return 0.0;
	}

	double get_z(size_t id) const
	{
		if (id < resn)
			return pro.at(id).z;
		else
			return 0.0;
	}

	double get_bfactor(size_t id) const
	{
		if (id < resn)
			return pro.at(id).bfactor;
		else
			return 0.0;
	}

	VectorXd get_rescoord(size_t id) const
	{
		if (id < resn)
			return rescoords.at(id);
		else
			return VectorXd();
	}

	size_t get_resatomn(size_t id) const
	{
		if (id < resn)
			return rescoords.at(id).size() / 3;
		else
			return 0;
	}

	size_t get_resn() const { return resn; }

	int get_contact(size_t i, size_t j) const
	{
		if (i < resn && j < resn)
			return contact_map(i, j);
		else
			return 0;
	}

	MatrixXi get_contact_map() const { return contact_map; }

	VectorXd get_procoord() const { return procoord; }

	VectorXd get_ligandcoord() const { return ligandcoord; }

	VectorXd get_dist2ligand() const { return dist2ligand; }

	MatrixXd get_distmat() const { return distmat; }

	ArrayXXd get_kmat() const { return kmat; }

	MatrixXd gen_hessian();

	MatrixXd gen_covariance(MatrixXd hessian);

	void gen_entropy(MatrixXd hessian, double &entropy);

	void show_contact_pairs();

private:
	double distance(size_t i, size_t j);

	size_t calc_zero_modes(VectorXd eigenvalues, std::vector<size_t> *zeromodes, std::vector<size_t> *nonzeromodes);

	size_t calc_zero_modes(VectorXd eigenvalues, VectorXd &zero2inf_eigenvalues);

	void read(std::string fpath);

	void gen_contact();

	void gen_coord();

	void gen_distmat();

	void gen_dist2ligand();

	std::map<size_t, ResInfo> pro;
	std::map<size_t, AtomInfoList> proatoms;
	std::map<std::string, AtomInfoList> ligand;
	std::map<std::string, AtomInfoList> excl;

	MatrixXi contact_map;
	std::vector<ContactPair> contact_pairs;
	ArrayXXd kmat;
	MatrixXd distmat;
	VectorXd procoord;
	std::map<size_t, VectorXd> rescoords;
	VectorXd ligandcoord;
	VectorXd dist2ligand;

	std::set<std::string> prores = { "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" };
	std::set<std::string> ligandres = {};
	std::set<std::string> exclres = { "HOH" };

	size_t resn = 0;
	size_t proatomn = 0;
	size_t ligandatomn = 0;
	size_t pairn = 0;

	double cutoff_intra = 0.0;
	double cutoff_inter = 0.0;
	double k_intra = 0.0;
	double k_inter = 0.0;

	bool gen_distmat_flag = false;
	bool gen_contact_flag = false;

	bool with_ligand_flag = false;
};
