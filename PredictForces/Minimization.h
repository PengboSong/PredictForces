#pragma once
#include "defines.h"
#include "method.h"

extern HandleMessage Log;

class Minimization
{
public:
	Minimization();
	Minimization(BGDpara paras, size_t N, VectorXd apo_procoord, ArrayXXd k_matrix);
	~Minimization();

	bool convergence() const { return convergeStatus; }

	bool complete() const { return completeFlag; }

	VectorXd equilibrium_coord() { return coord; }

	VectorXd pocket_force() { return pocketforce; }

	void reset();

	void single_pocket(PocketList pocket_members, VectorXd fix_coord);

	void single_pocket(VectorXd pocket_force);

	void single_pocket(PocketList pocket_members, VectorXd pocket_force, VectorXd fix_coord);

	void multiple_pocket(PocketList fixed_pocket_members, VectorXd fixed_pocket_coord, VectorXd fixed_pocket_force, PocketList rotate_pocket_members, VectorXd rotate_pocket_coord, VectorXd rotate_pocket_force, Vector3d & vcenter, Vector3d & degrees);

	FreeEnergy energy();

	void print_gradient() { print_xyzvec(dR, "Gradient on residue"); }

	void print_coord() { print_xyzvec(coord, "XYZ coordinates of residue"); }

	void print_pocketforce() { print_xyzvec(pocketforce, "Force on residue"); }
	
private:
	void update_coords_derivate();

	Vector3d euler_degrees_derivate(Vector3d r, Vector3d dvdr, Vector3d degrees);

	VectorXd degrees_derivate(VectorXd domain_coord, Vector3d vcenter, Vector3d degrees);

	void print_xyzvec(VectorXd vec, std::string prefix = std::string());

	bool convergeStatus = false;
	bool completeFlag = false;
	bool downwardStep = false;
	BGDpara parameters;
	double realstep = 0.0;

	size_t noden = 0, coordn = 0;

	VectorXd coord, coord_apo, coord_ref;
	VectorXd pocketforce;
	MatrixXd distmat, distmat_apo, distmat_ref;
	MatrixXd xyzdiffmat;
	ArrayXXd kmat;
	MatrixXd coeff;

	double potential = 0.0, prev_potential = 0.0, ref_potential = 0.0;
	double err = 0.0;

	VectorXd dR;
};

