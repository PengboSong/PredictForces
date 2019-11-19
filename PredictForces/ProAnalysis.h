#pragma once
#include <fstream>
#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include "defines.h"
#include "method.h"
#include "Pro.h"
#include "Minimization.h"

extern HandleMessage Log;

class ProAnalysis
{
public:
	ProAnalysis();
	ProAnalysis(Pro apo, Pro binding, Pro allostery, Pro complex);
	~ProAnalysis();

	void interactive_pocket(Pockets m);

	void interactive();

	MatrixXd get_hessian() const { return hessian; }

	Matrix3d get_hessian(size_t i, size_t j);

	double get_hessian_s(size_t si, size_t sj);

	MatrixXd get_covariance() { return covariance; }

	Matrix3d get_covariance(size_t i, size_t j);

	double get_covariance_s(size_t si, size_t sj);

	void show_LFmethod_detail();

	void show_pocket(Pockets m) { show_pocket(pocket_members(m)); }

	void show_pocket_force(Pockets m) { show_pocket_force(pocket(m).access_force, pocket(m).members, pocket(m).force); }

	void show_pro_pocket_force(Pockets m) { show_pro_pocket_force(pocket_members(m), pro(m).force); }

	void show_pro_all_force(Pockets m) { show_pro_all_force(pro(m).force); }

	void test_pocket(Pockets m)
	{
		test_pocket(
			pocket(m).access_force,
			pro(m).preprocess,
			pocket(m).members,
			pro(m).displacement,
			pocket(m).force,
			pro(m).fitcoord
		);
		calc_model_correlation(
			pocket(m).access_force,
			pocket(m).force,
			pro(m).displacement
		);
	}

	bool in_pocket(Pockets m, size_t id) { return in_pocket(pocket_members(m), id); }

	void add_to_pocket(Pockets m, size_t id) { add_to_pocket(pocket_members(m), id); }

	void remove_from_pocket(Pockets m, size_t id) { remove_from_pocket(pocket_members(m), id); }

	void gen_pocket(Pockets m, double cutoff)
	{
		gen_pocket(
			proinfos.at(m).withligand,
			pockets.at(m).members,
			cutoff,
			proinfos.at(m).dist2ligand
		);
	}

	double calc_apo_pro_rmsd(Pockets m)
	{
		if (pro(m).preprocess)
			calc_model_rmsd(
				pocket(m).access_force,
				pocket(m).force,
				pro(m).fitcoord
			);
	}

	void gen_pocket_force(Pockets m);

	void gen_free_energy();

	void write_matrix(MatrixXd mat, std::string writepath);

	void write_matrix_binary(MatrixXd mat, std::string writepath);

	void read_matrix_binary(MatrixXd & mat, std::string fpath);

private:
	std::string write_path(std::string fname);

	PocketInfo & pocket(Pockets m) { return pockets.at(m); }

	PocketList & pocket_members(Pockets m) { return pocket(m).members; }

	ProInfo & pro(Pockets m) { return proinfos.at(m); }

	void init_container();

	void preprocess(Pockets m);

	void switch_LFmethod(VectorXd &coeff, MatrixXd X, VectorXd Y);

	void energy_minimization(EMType em, Pockets m);

	double calc_model_rmsd(bool access, VectorXd pocket_force, VectorXd refcoord);

	double calc_model_correlation(bool access, VectorXd pocket_force, VectorXd displacement);

	double calc_correlation(VectorXd displacement, VectorXd new_displacement);

	void show_pocket(PocketList pocket_members);

	void test_pocket(bool access, bool preprocess, PocketList pocket_members, VectorXd displacement, VectorXd pocket_force, VectorXd refcoord);

	void show_pocket_force(bool access, PocketList pocket_members, VectorXd pocket_force);

	void show_pro_pocket_force(PocketList pocket_members, VectorXd pro_force);

	void show_pro_all_force(VectorXd pro_force);

	bool in_pocket(PocketList pocket_members, size_t id);

	void add_to_pocket(PocketList & pocket_members, size_t id);

	void remove_from_pocket(PocketList & pocket_members, size_t id);

	void gen_pocket(bool has_ligand, PocketList & pocket_members, double cutoff, VectorXd dist2ligand);

	void gen_pocket_force(bool & access, VectorXd & pocket_force, PocketList pocket_members, VectorXd pro_force, VectorXd displacement);

	void gen_pocket_force(bool & access, VectorXd & pocket_force, VectorXd fixed_force, PocketList pocket_members, PocketList fixed_pocket, VectorXd pro_force, VectorXd displacement);

	// Pack functions for convenience
	void check_EM_diff(Pockets m);

	void print_energy_results();

	Pro ProE, ProS, ProA, ProAS;

	PocketContainer pockets;

	ProInfoContainer proinfos;

	MatrixXd hessian;
	MatrixXd covariance;

	VectorXd apo_procoord;

	MatrixXd mprocoord;
	
	double ddG = 0.0;
	double ddG_predict = 0.0;
	double ddG_deduce = 0.0;
	double ddG_deduce1 = 0.0;
	double ddG_deduce2 = 0.0;

	// Multiple linear fitting method
	LFmethods LFMETHOD = BatchGradientDescent;
	Minimization MinPocket;
	BGDpara BGDparameters = { 1e-2, 1e-4, 1.2, 1.5, 1000000 };

	double cutoff_intra = 9.0;
	double k_intra = 10.0;
	double pocket_cutoff = 4.5;

	std::string workdirPath;
};
