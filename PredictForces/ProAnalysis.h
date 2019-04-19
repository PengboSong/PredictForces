#pragma once
#include <algorithm>
#include <list>
#include <utility>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include "handle_io.h"
#include "Pro.h"

using namespace std;
using namespace Eigen;
using boost::format;
using boost::lexical_cast;
using boost::algorithm::trim;

constexpr double PI = 3.1415926535897932;

class ProAnalysis
{
public:
	ProAnalysis();
	ProAnalysis(Pro apo, Pro binding);
	ProAnalysis(Pro apo, Pro binding, Pro allostery, Pro complex);
	~ProAnalysis();

	void interactive_pocket(unsigned int mode);
	void interactive();

	MatrixXd get_hessian();
	Matrix3d get_hessian(size_t i, size_t j);
	double get_hessian_s(size_t si, size_t sj);
	void write_hessian(string writepath);

	MatrixXd get_covariance();
	Matrix3d get_covariance(size_t i, size_t j);
	double get_covariance_s(size_t si, size_t sj);
	void write_covariance(string writepath);

	void set_learning_step(double step)
	{
		if (step > 0)
			LEARNING_STEP = step;
	}
	void set_convergence(double limit)
	{
		if (limit > 0)
			CONVERGENCE = limit;
	}
	void set_iteration_times(size_t N)
	{
		ITERATION_TIMES = N;
	}
	void set_random_times(size_t N)
	{
		RANDOM_TIMES = N;
	}

	void show_LFmethod_detail();

	void set_LFmethod(unsigned int mode);

	void choose_LFmethod();

	list<size_t> get_pocketS() {
		return pocketS;
	}
	list<size_t> get_pocketA() {
		return pocketA;
	}
	list<size_t> get_pocketAS() {
		return pocketAS;
	}

	void show_pocketS() {
		show_pocket(pocketS);
	}
	void show_pocketA() {
		show_pocket(pocketA);
	}
	void show_pocketAS() {
		show_pocket(pocketAS);
	}

	void show_pocketS_force() {
		show_pocket_force(has_pocketS_force_flag, pocketS, pocketS_force);
	}
	void show_pocketA_force() {
		show_pocket_force(has_pocketA_force_flag, pocketA, pocketA_force);
	}
	void show_pocketAS_force() {
		show_pocket_force(has_pocketAS_force_flag, pocketAS, pocketAS_force);
	}

	void show_pro_pocketS_force() {
		show_pro_pocket_force(pocketS, ES_force);
	}
	void show_pro_pocketA_force() {
		show_pro_pocket_force(pocketA, EA_force);
	}
	void show_pro_pocketAS_force() {
		show_pro_pocket_force(pocketAS, EAS_force);
	}

	void show_proS_all_force() {
		show_pro_all_force(ES_force);
	}
	void show_proA_all_force() {
		show_pro_all_force(EA_force);
	}
	void show_proAS_all_force() {
		show_pro_all_force(EAS_force);
	}

	void test_pocketS() {
		test_pocket(has_pocketS_force_flag, ES_info, pocketS, pocketS_force, S_fitprocoord);
	}
	void test_pocketA() {
		test_pocket(has_pocketA_force_flag, EA_info, pocketA, pocketA_force, A_fitprocoord);
	}
	void test_pocketAS() {
		test_pocket(has_pocketAS_force_flag, EAS_info, pocketAS, pocketAS_force, AS_fitprocoord);
	}

	bool in_pocketS(size_t id) {
		return in_pocket(pocketS, id);
	}
	bool in_pocketA(size_t id) {
		return in_pocket(pocketA, id);
	}
	bool in_pocketAS(size_t id) {
		return in_pocket(pocketAS, id);
	}

	void add_to_pocketS(size_t id) {
		add_to_pocket(pocketS, id);
	}
	void add_to_pocketA(size_t id) {
		add_to_pocket(pocketA, id);
	}
	void add_to_pocketAS(size_t id) {
		add_to_pocket(pocketAS, id);
	}

	void remove_from_pocketS(size_t id) {
		remove_from_pocket(pocketS, id);
	}
	void remove_from_pocketA(size_t id) {
		remove_from_pocket(pocketA, id);
	}
	void remove_from_pocketAS(size_t id) {
		remove_from_pocket(pocketAS, id);
	}

	void gen_pocketS(double cutoff) {
		if (!ProS.empty())
			gen_pocket(ProS.has_ligand(), pocketS, cutoff, S_dist2ligand);
	}
	void gen_pocketA(double cutoff) {
		if (!ProA.empty())
			gen_pocket(ProA.has_ligand(), pocketA, cutoff, A_dist2ligand);
	}
	void gen_pocketAS(double cutoff) {
		if (!ProAS.empty())
			gen_pocket(ProAS.has_ligand(), pocketAS, cutoff, AS_dist2ligand);
	}

	double calc_ES_rmsd() {
		if (ES_info)
			calc_model_rmsd(has_pocketS_force_flag, pocketS_force, S_fitprocoord);
	}
	double calc_EA_rmsd() {
		if (EA_info)
			calc_model_rmsd(has_pocketA_force_flag, pocketA_force, A_fitprocoord);
	}
	double calc_EAS_rmsd() {
		if (EAS_info)
			calc_model_rmsd(has_pocketAS_force_flag, pocketAS_force, AS_fitprocoord);
	}

	void gen_pocketS_force() {
		if (ES_info)
			gen_pocket_force(has_pocketS_force_flag, pocketS_force, pocketS, ES_force, ES_displacement);
	}
	void gen_pocketA_force() {
		if (EA_info)
			gen_pocket_force(has_pocketA_force_flag, pocketA_force, pocketA, EA_force, EA_displacement);
	}
	void gen_pocketAS_force() {
		if (EAS_info)
		{
			if (has_pocketS_force_flag)
				gen_pocket_force(has_pocketAS_force_flag, pocketAS_force, pocketS_force, pocketAS, pocketS, EAS_force, EAS_displacement);
			else if (has_pocketA_force_flag)
				gen_pocket_force(has_pocketAS_force_flag, pocketAS_force, pocketA_force, pocketAS, pocketA, EAS_force, EAS_displacement);
			else
				gen_pocket_force(has_pocketAS_force_flag, pocketAS_force, pocketAS, EAS_force, EAS_displacement);
		}
	}

	void gen_free_energy();

private:
	void switch_LFmethod(VectorXd &coeff, MatrixXd X, MatrixXd Y);

	double calc_model_rmsd(bool flag, VectorXd pocket_force, VectorXd refcoord);

	void show_pocket(list<size_t> pocket);

	void test_pocket(bool flag, bool info, list<size_t> pocket, VectorXd pocket_force, VectorXd refcoord);

	void show_pocket_force(bool flag, list<size_t> pocket, VectorXd pocket_force);

	void show_pro_pocket_force(list<size_t> pocket, VectorXd pro_force);

	void show_pro_all_force(VectorXd pro_force);

	bool in_pocket(list<size_t> pocket, size_t id);

	void add_to_pocket(list<size_t> & pocket, size_t id);

	void remove_from_pocket(list<size_t> & pocket, size_t id);

	void gen_pocket(bool has_ligand, list<size_t> & pocket, double cutoff, VectorXd dist2ligand);

	void gen_pocket_force(bool & flag, VectorXd & pocket_force, list<size_t> pocket, VectorXd pro_force, VectorXd displacement);

	void gen_pocket_force(bool & flag, VectorXd & pocket_force, VectorXd fixed_force, list<size_t> pocket, list<size_t> fixed_pocket, VectorXd pro_force, VectorXd displacement);

	void calc_energy_known(bool flag, double & proenergy, double & pocketenergy, double & totenergy, list<size_t> pocket, VectorXd pro_force, MatrixXd distmat);

	void calc_energy_unknown(bool flag, double & proenergy, double & pocketenergy, double & totenergy, VectorXd pocket_force);

	void print_energy_results();

	Pro ProE;
	Pro ProS;
	Pro ProA;
	Pro ProAS;

	VectorXd S_fitprocoord;
	VectorXd S_dist2ligand;
	VectorXd pocketS_force;
	double S_proenergy = 0.0;
	double S_pocketenergy = 0.0;
	double S_energy = 0.0;
	double S_predict_proenergy = 0.0;
	double S_predict_pocketenergy = 0.0;
	double S_predict_energy = 0.0;

	VectorXd A_fitprocoord;
	VectorXd A_dist2ligand;
	VectorXd pocketA_force;
	double A_proenergy = 0.0;
	double A_pocketenergy = 0.0;
	double A_energy = 0.0;
	double A_predict_proenergy = 0.0;
	double A_predict_pocketenergy = 0.0;
	double A_predict_energy = 0.0;

	VectorXd AS_fitprocoord;
	VectorXd AS_dist2ligand;
	VectorXd pocketAS_force;
	double AS_proenergy = 0.0;
	double AS_pocketenergy = 0.0;
	double AS_energy = 0.0;
	double AS_predict_proenergy = 0.0;
	double AS_predict_pocketenergy = 0.0;
	double AS_predict_energy = 0.0;

	MatrixXd hessian;
	MatrixXd covariance;

	double ddG = 0.0;
	double ddG_predict = 0.0;

	bool ES_info = false;
	VectorXd ES_displacement;
	VectorXd ES_force;
	double ES_average_force;
	double ES_rmsd;

	bool EA_info = false;
	VectorXd EA_displacement;
	VectorXd EA_force;
	double EA_average_force;
	double EA_rmsd;

	bool EAS_info = false;
	VectorXd EAS_displacement;
	VectorXd EAS_force;
	double EAS_average_force;
	double EAS_rmsd;

	list<size_t> pocketS;
	list<size_t> pocketA;
	list<size_t> pocketAS;

	// Status
	bool has_pocketS_force_flag = false;
	bool has_pocketA_force_flag = false;
	bool has_pocketAS_force_flag = false;

	// Matrix formats
	IOFormat CleanFmt = IOFormat(4, 0, ", ", "\n", "[", "]");

	// Multiple linear fitting method
	unsigned int LFmethod_mode = 0;
	map<unsigned int, string> LFmethods = { {0, "Normal Equation"}, {1, "Batch Gradient Descent"} };

	// BGD parameters
	double LEARNING_STEP = 1e-6;
	double CONVERGENCE = 1e-2;
	size_t ITERATION_TIMES = 10000;
	size_t RANDOM_TIMES = 1000;
};
