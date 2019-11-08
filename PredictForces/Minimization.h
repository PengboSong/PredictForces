#pragma once
#include "defines.h"
#include "method.h"

class Minimization
{
public:
	Minimization();
	Minimization(BGDpara paras, size_t N, VectorXd apo_procoord, ArrayXXd k_matrix);
	~Minimization();

	bool convergence() const { return convergeStatus; }

	void clear();

	void single_pocket(PocketList pocket_members, VectorXd & pocket_force, VectorXd fix_coord, VectorXd & equilibrium_coord);

	void single_pocket(VectorXd & pocket_force, VectorXd & equilibrium_coord);

	void single_pocket_with_force(PocketList pocket_members, VectorXd & pocket_force, VectorXd fix_coord, VectorXd & equilibrium_coord);

	void multiple_pocket(PocketList fixed_pocket_members, VectorXd fixed_pocket_coord, PocketList rotate_pocket_members, VectorXd rotate_pocket_coord, VectorXd & equilibrium_coord, Vector3d & vcenter, Vector3d & degrees);
	
private:
	void update_coords_derivate();

	Vector3d euler_degrees_derivate(Vector3d r, Vector3d dvdr, Vector3d degrees);

	VectorXd degrees_derivate(VectorXd domain_coord, Vector3d vcenter, Vector3d degrees);

	bool convergeStatus = false;
	bool downwardStep = false;
	BGDpara parameters;
	double realstep = 0.0;

	size_t noden = 0, coordn = 0;

	VectorXd coord, coord_apo;
	MatrixXd distmat, distmat_apo, distmat_ref;
	MatrixXd xyzdiffmat;
	ArrayXXd kmat;

	MatrixXd coeff;

	double potential = 0.0, prev_potential = 0.0, ref_potential = 0.0;
	double err = 0.0;

	VectorXd dR;
};

