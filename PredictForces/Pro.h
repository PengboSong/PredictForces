#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>

#include <Eigen/Dense>

#include "read_pdb.h"

class Pro
{
public:
	Pro();
	Pro(std::string fpath, bool has_ligand_flag, std::set<std::string> exclude = { }, double k = 10.0, double cutoff = 1.0);
	~Pro();

	bool has_ligand();

	bool has_res(size_t id);
	std::string get_resname(size_t id);
	std::string get_chain(size_t id);
	size_t get_resid(size_t id);
	double get_x(size_t id);
	double get_y(size_t id);
	double get_z(size_t id);
	double get_bfactor(size_t id);

	Eigen::VectorXd get_rescoord(size_t id);
	size_t get_resatomn(size_t id);

	size_t get_resn();

	int get_contact(size_t i, size_t j);

	Eigen::VectorXd get_procoord();
	Eigen::VectorXd get_ligandcoord();

	Eigen::VectorXd get_dist2ligand();

	Eigen::MatrixXd get_distmat();
	Eigen::ArrayXXd get_kmat();

	Eigen::MatrixXd get_hessian();
	Eigen::Matrix3d get_hessian(size_t i, size_t j);
	double get_hessian_s(size_t si, size_t sj);

	void write_hessian(std::string writepath);
	void write_covariance(std::string writepath);

	bool empty();

private:
	double distance(size_t i, size_t j);
	double diff_x(size_t i, size_t j);
	double diff_y(size_t i, size_t j);
	double diff_z(size_t i, size_t j);

	double distance(pair ij);
	double diff_x(pair ij);
	double diff_y(pair ij);
	double diff_z(pair ij);

	void read(std::string fpath);

	void gen_contact();
	
	void gen_coord();

	void gen_hessian();

	void gen_covariance();

	void gen_distmat();

	void gen_dist2ligand();

	std::map<size_t, ResInfo> pro;
	std::map<size_t, std::vector<AtomInfo>> proatoms;
	std::map<std::string, std::vector<AtomInfo>> ligand;
	std::map<std::string, std::vector<AtomInfo>> excl;

	Eigen::MatrixXi contact_map;
	std::vector<pair> contact_pairs;
	Eigen::ArrayXXd kmat;
	Eigen::MatrixXd distmat;
	Eigen::MatrixXd hessian;
	Eigen::MatrixXd covariance;
	Eigen::VectorXd procoord;
	std::map<size_t, Eigen::VectorXd> rescoords;
	Eigen::VectorXd ligandcoord;
	Eigen::VectorXd dist2ligand;

	std::set<std::string> prores = { "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL" };
	std::set<std::string> ligandres = {};
	std::set<std::string> exclres = { "HOH" };

	Eigen::IOFormat CleanFmt = Eigen::IOFormat(4, 0, ", ", "\n", "[", "]");

	size_t resn = 0;
	size_t proatomn = 0;
	size_t ligandatomn = 0;
	size_t pairn = 0;
	double cutoff_intra = 10.0;
	double cutoff_inter = 10.0;
	double k_intra = 1.0;
	double k_inter = 1.0;
	double k_default = 1.0;

	bool gen_distmat_flag = false;
	bool gen_contact_flag = false;
	bool gen_hessian_flag = false;
	bool gen_covariance_flag = false;

	bool with_ligand_flag = false;
};
