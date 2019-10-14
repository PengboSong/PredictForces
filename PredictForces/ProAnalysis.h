#pragma once
#include <algorithm>
#include <fstream>
#include <list>
#include <utility>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include "handle_io.h"
#include "Pro.h"
#include "method.h"

using namespace Eigen;

namespace filesys = boost::filesystem;

constexpr double PI = 3.1415926535897932;

enum Pockets : uint8_t {
	POCKETS,	  // Binding Pocket
	POCKETA,      // Allostery Pocket
	POCKETAS	  // Complex Pocket
};

enum LFmethods : uint8_t {
	BatchGradientDescent,
	NormalEquation
};

typedef std::list<size_t> PocketList;

typedef struct {
	bool access_force = false;
	PocketList members;
	VectorXd force;
} PocketInfo;


typedef struct {
	double pro = 0.0;
	double pocket = 0.0;
	double total = 0.0;
} FreeEnergy;

typedef struct {
	bool empty = true;
	bool withligand = false;
	VectorXd procoord;
	VectorXd fitprocoord;
	VectorXd dist2ligand;
	FreeEnergy G;
} ProInfo;

typedef struct {
	bool preprocess = false;
	double rmsd = 0.0;
	double meanforce = 0.0;
	double pearson = 0.0;
	VectorXd displacement;
	VectorXd force;
	VectorXd equilibrium_coord;
	MatrixXd distdiff;
} ApoProInfo;

typedef std::map<Pockets, PocketInfo> PocketContainer;
typedef std::map<Pockets, ProInfo> ProInfoContainer;
typedef std::map<Pockets, ApoProInfo> ApoProInfoContainer;

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

	MatrixXd gen_full_distmat(VectorXd coord);
	MatrixXd gen_all_distmat(VectorXd coord);
	MatrixXd gen_small_distmat(VectorXd coord);

	void show_LFmethod_detail();

	void show_pocket(Pockets m) { show_pocket(pocket_members(m)); }

	void show_pocket_force(Pockets m) { show_pocket_force(pocket(m).access_force, pocket(m).members, pocket(m).force); }

	void show_pro_pocket_force(Pockets m) { show_pro_pocket_force(pocket_members(m), apo_pro(m).force); }

	void show_pro_all_force(Pockets m) { show_pro_all_force(apo_pro(m).force); }

	void test_pocket(Pockets m)
	{
		test_pocket(
			pocket(m).access_force,
			apo_pro(m).preprocess,
			pocket(m).members,
			apo_pro(m).displacement,
			pocket(m).force,
			pro(m).fitprocoord
		);
		calc_model_correlation(
			pocket(m).access_force,
			pocket(m).force,
			apo_pro(m).displacement
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
		if (apo_pro(m).preprocess)
			calc_model_rmsd(
				pocket(m).access_force,
				pocket(m).force,
				pro(m).fitprocoord
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

	ApoProInfo & apo_pro(Pockets m) { return apo_proinfos.at(m); }

	void init_container();

	void preprocess(Pockets m);

	void switch_LFmethod(VectorXd &coeff, MatrixXd X, VectorXd Y);

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

	void calc_energy_known(bool access, FreeEnergy & energy, PocketList pocket_members, VectorXd pro_force, MatrixXd distmat, VectorXd displacement);

	void calc_energy_unknown(bool access, FreeEnergy & energy, VectorXd pocket_force, VectorXd equilibrium_coord);

	void minimization(bool & flag, PocketList pocket_members, VectorXd fix_procoord, VectorXd & equilibrium_coord, VectorXd & pocket_force, BGDpara paras);

	void minimization(bool & flag, VectorXd & equilibrium_coord, VectorXd & pocket_force, BGDpara paras);

	void optimize_pocket_structure(Vector3d & vcenter, Vector3d & degrees, Matrix3Xd coord, MatrixXd distmat_0, MatrixXd distmat, MatrixXi contactmap, double k);

	// Pack functions for convenience
	void minimization_calc_energy(bool & flag, Pockets m);

	void equilibrium_coord_rmsd(Pockets m);

	void print_energy_results();

	Pro ProE, ProS, ProA, ProAS;

	double ES_pearson;
	double EA_pearson;
	double EAS_pearson;

	PocketContainer pockets;

	ProInfoContainer proinfos;

	ApoProInfoContainer apo_proinfos;

	MatrixXd hessian;
	MatrixXd covariance;

	VectorXd apo_procoord;

	MatrixXd mprocoord;
	
	double ddG = 0.0;
	double ddG_predict = 0.0;
	double ddG_deduce = 0.0;
	double ddG_deduce1 = 0.0;
	double ddG_deduce2 = 0.0;

	// Matrix formats
	IOFormat CleanFmt = IOFormat(4, 0, ", ", "\n", "[", "]");

	// Multiple linear fitting method
	LFmethods LFMETHOD = BatchGradientDescent;

	// BGD parameters
	BGDpara BGDparameters = { 1e-2, 1e-4, 1000000 };

	double cutoff_intra = 9.0;
	double k_intra = 10.0;
	double pocket_cutoff = 4.5;

	std::string workdir_path;
};
