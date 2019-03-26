#pragma once
#include <iomanip>
#include <algorithm>
#include <list>

#include "Pro.h"

constexpr double PI = 3.1415926535897932;
constexpr double Navo = 6.02214076e23;
constexpr double kB = 1.380649e-23;
constexpr double Temp = 298.15;

class ProAnalysis
{
public:
	ProAnalysis();
	ProAnalysis(Pro apo, Pro binding);
	ProAnalysis(Pro apo, Pro binding, Pro allostery, Pro complex);
	~ProAnalysis();

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

	void test_pocketS() {
		test_pocket(ES_info, pocketS, pocketS_force, ES_force, ProS.get_procoord());
	}
	void test_pocketA() {
		test_pocket(EA_info, pocketA, pocketA_force, EA_force, ProA.get_procoord());
	}
	void test_pocketAS() {
		test_pocket(EAS_info, pocketAS, pocketAS_force, EAS_force, ProAS.get_procoord());
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
		gen_pocket(pocketS, cutoff, S_dist2ligand);
	}
	void gen_pocketA(double cutoff) {
		gen_pocket(pocketA, cutoff, A_dist2ligand);
	}
	void gen_pocketAS(double cutoff) {
		gen_pocket(pocketAS, cutoff, AS_dist2ligand);
	}

	double calc_ES_rmsd() {
		calc_model_rmsd(pocketS, pocketS_force, ES_force, ProS.get_procoord());
	}
	double calc_EA_rmsd() {
		calc_model_rmsd(pocketA, pocketA_force, EA_force, ProA.get_procoord());
	}
	double calc_EAS_rmsd() {
		calc_model_rmsd(pocketAS, pocketAS_force, EAS_force, ProAS.get_procoord());
	}

	void gen_pocketS_force() {
		gen_pocket_force(pocketS_force, pocketS, ES_force);
	}
	void gen_pocketA_force() {
		gen_pocket_force(pocketA_force, pocketA, EA_force);
	}
	void gen_pocketAS_force() {
		gen_pocket_force(pocketAS_force, pocketAS, EAS_force);
	}

	void gen_free_energy();
	
private:


	double calc_model_rmsd(std::list<size_t> pocket, Eigen::VectorXd pocket_force, Eigen::VectorXd pro_force, Eigen::VectorXd refcoord);
	
	void show_pocket(std::list<size_t> pocket);

	void test_pocket(bool info, std::list<size_t> pocket, Eigen::VectorXd pocket_force, Eigen::VectorXd pro_force, Eigen::VectorXd refcoord);

	void show_pocket_force(std::list<size_t> pocket, Eigen::VectorXd pocket_force);

	bool in_pocket(std::list<size_t> pocket, size_t id);

	void add_to_pocket(std::list<size_t> pocket, size_t id);

	void remove_from_pocket(std::list<size_t> pocket, size_t id);

	void gen_pocket(std::list<size_t> &pocket, double cutoff, Eigen::VectorXd dist2ligand);

	void gen_pocket_force(Eigen::VectorXd & pocket_force, std::list<size_t> pocket, Eigen::VectorXd pro_force);

	void calc_energy_known(double & proenergy, double & pocketenergy, double & totenergy, Eigen::VectorXd pocket_force, Eigen::MatrixXd distmat);

	void calc_energy_unknown(double & proenergy, double & pocketenergy, double & totenergy, Eigen::VectorXd pocket_force);

	void debug_energy_unknown(Eigen::VectorXd pocket_force);
		
	void print_energy_results();

	Pro ProE;
	Pro ProS;
	Pro ProA;
	Pro ProAS;

	Eigen::VectorXd S_dist2ligand;
	Eigen::VectorXd pocketS_force;
	double S_mincutoff = 0.0;
	double S_proenergy = 0.0;
	double S_pocketenergy = 0.0;
	double S_energy = 0.0;

	Eigen::VectorXd A_dist2ligand;
	Eigen::VectorXd pocketA_force;
	double A_mincutoff = 0.0;
	double A_proenergy = 0.0;
	double A_pocketenergy = 0.0;
	double A_energy = 0.0;

	Eigen::VectorXd AS_dist2ligand;
	Eigen::VectorXd pocketAS_force;
	double AS_mincutoff = 0.0;
	double AS_proenergy = 0.0;
	double AS_pocketenergy = 0.0;
	double AS_energy = 0.0;

	Eigen::MatrixXd H;
	Eigen::MatrixXd G;

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

	size_t MAX_SCAN_RESN = 20;
	double CUTOFF = 10.0;
};

