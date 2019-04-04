#pragma once
#include <iomanip>
#include <algorithm>
#include <list>

#include "Pro.h"

constexpr double PI = 3.1415926535897932;

class ProAnalysis
{
public:
	ProAnalysis();
	ProAnalysis(Pro apo, Pro binding);
	ProAnalysis(Pro apo, Pro binding, Pro allostery, Pro complex);
	~ProAnalysis();

	Eigen::MatrixXd get_hessian();
	Eigen::Matrix3d get_hessian(size_t i, size_t j);
	double get_hessian_s(size_t si, size_t sj);
	void write_hessian(std::string writepath);

	Eigen::MatrixXd get_covariance();
	Eigen::Matrix3d get_covariance(size_t i, size_t j);
	double get_covariance_s(size_t si, size_t sj);
	void write_covariance(std::string writepath);

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
	void set_iteration_times(long N)
	{
		if (N > 0)
			ITERATION_TIMES = size_t(N);
	}

	std::list<size_t> get_pocketS() {
		return pocketS;
	}
	std::list<size_t> get_pocketA() {
		return pocketA;
	}
	std::list<size_t> get_pocketAS() {
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
		show_pocket_force(pocketS, pocketS_force);
	}
	void show_pocketA_force() {
		show_pocket_force(pocketA, pocketA_force);
	}
	void show_pocketAS_force() {
		show_pocket_force(pocketAS, pocketAS_force);
	}

	void test_pocketS() {
		if (ES_info)
			test_pocket(ES_info, pocketS, pocketS_force, ES_force, S_fitprocoord, ES_displacement);
	}
	void test_pocketA() {
		if (EA_info)
			test_pocket(EA_info, pocketA, pocketA_force, EA_force, A_fitprocoord, EA_displacement);
	}
	void test_pocketAS() {
		if (EAS_info)
			test_pocket(EAS_info, pocketAS, pocketAS_force, EAS_force, AS_fitprocoord, EAS_displacement);
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
			calc_model_rmsd(pocketS, pocketS_force, ES_force, S_fitprocoord, ES_displacement);
	}
	double calc_EA_rmsd() {
		if (EA_info)
			calc_model_rmsd(pocketA, pocketA_force, EA_force, A_fitprocoord, EA_displacement);
	}
	double calc_EAS_rmsd() {
		if (EAS_info)
			calc_model_rmsd(pocketAS, pocketAS_force, EAS_force, AS_fitprocoord, EAS_displacement);
	}

	void gen_pocketS_force() {
		if (ES_info)
			gen_pocket_force(pocketS_force, pocketS, ES_force, ES_displacement);
	}
	void gen_pocketA_force() {
		if (EA_info)
			gen_pocket_force(pocketA_force, pocketA, EA_force, EA_displacement);
	}
	void gen_pocketAS_force() {
		if (EAS_info)
			gen_pocket_force(pocketAS_force, pocketAS, EAS_force, EAS_displacement);
	}

	void gen_free_energy();
	
private:
	double calc_model_rmsd(std::list<size_t> pocket, Eigen::VectorXd pocket_force, Eigen::VectorXd pro_force, Eigen::VectorXd refcoord, Eigen::VectorXd displacenment);
	
	void show_pocket(std::list<size_t> pocket);

	void test_pocket(bool info, std::list<size_t> pocket, Eigen::VectorXd pocket_force, Eigen::VectorXd pro_force, Eigen::VectorXd refcoord, Eigen::VectorXd displacenment);

	void show_pocket_force(std::list<size_t> pocket, Eigen::VectorXd pocket_force);

	bool in_pocket(std::list<size_t> pocket, size_t id);

	void add_to_pocket(std::list<size_t> pocket, size_t id);

	void remove_from_pocket(std::list<size_t> pocket, size_t id);

	void gen_pocket(bool has_ligand, std::list<size_t> &pocket, double cutoff, Eigen::VectorXd dist2ligand);

	void gen_pocket_force(Eigen::VectorXd & pocket_force, std::list<size_t> pocket, Eigen::VectorXd pro_force, Eigen::VectorXd displacenment);
	
	void calc_energy_known(double & proenergy, double & pocketenergy, double & totenergy, Eigen::VectorXd pocket_force, Eigen::MatrixXd distmat);

	void calc_energy_unknown(double & proenergy, double & pocketenergy, double & totenergy, Eigen::VectorXd pocket_force);

	void debug_energy_unknown(Eigen::VectorXd pocket_force);
		
	void print_energy_results();

	Pro ProE;
	Pro ProS;
	Pro ProA;
	Pro ProAS;

	Eigen::VectorXd S_fitprocoord;
	Eigen::VectorXd S_dist2ligand;
	Eigen::VectorXd pocketS_force;
	double S_proenergy = 0.0;
	double S_pocketenergy = 0.0;
	double S_energy = 0.0;

	Eigen::VectorXd A_fitprocoord;
	Eigen::VectorXd A_dist2ligand;
	Eigen::VectorXd pocketA_force;
	double A_proenergy = 0.0;
	double A_pocketenergy = 0.0;
	double A_energy = 0.0;

	Eigen::VectorXd AS_fitprocoord;
	Eigen::VectorXd AS_dist2ligand;
	Eigen::VectorXd pocketAS_force;
	double AS_proenergy = 0.0;
	double AS_pocketenergy = 0.0;
	double AS_energy = 0.0;

	Eigen::MatrixXd hessian;
	Eigen::MatrixXd covariance;

	double ddG = 0.0;

	bool ES_info = false;
	Eigen::VectorXd ES_displacement;
	Eigen::VectorXd ES_force;
	double ES_average_force;
	double ES_rmsd;

	bool EA_info = false;
	Eigen::VectorXd EA_displacement;
	Eigen::VectorXd EA_force;
	double EA_average_force;
	double EA_rmsd;

	bool EAS_info = false;
	Eigen::VectorXd EAS_displacement;
	Eigen::VectorXd EAS_force;
	double EAS_average_force;
	double EAS_rmsd;
	
	std::list<size_t> pocketS;
	std::list<size_t> pocketA;
	std::list<size_t> pocketAS;

	// Matrix formats
	Eigen::IOFormat CleanFmt = Eigen::IOFormat(4, 0, ", ", "\n", "[", "]");

	// BGD parameters
	double LEARNING_STEP = 1e-6;
	double CONVERGENCE = 1e-2;
	size_t ITERATION_TIMES = 10000;
};

