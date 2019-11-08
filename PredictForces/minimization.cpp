#include "ProAnalysis.h"

void ProAnalysis::minimization(bool & flag, PocketList pocket_members, VectorXd fix_coord, VectorXd & equilibrium_coord, VectorXd & pocket_force, BGDpara paras)
{
	double step = paras.learning_rate;

	// Set variables that would be often used
	VectorXd initial_coord = ProE.get_procoord();
	size_t veclen = initial_coord.size(); // 3N
	size_t resn = veclen / 3; // N
	size_t pocketlen = pocket_members.size() * 3; // L
	ArrayXXd kmat = ProE.get_kmat();
	MatrixXd dist_apo = gen_distmat(Dist, initial_coord);

	equilibrium_coord = initial_coord;
	modify_pocket_coord(equilibrium_coord, fix_coord, pocket_members);

	// Initialize variables
	double prev_potential = 0.0, potential = 0.0, ref_potential = 0.0, err = 0.0;
	// Pocket Force - L
	pocket_force = VectorXd::Zero(veclen);
	// Distance Matrix: Dist - N x N
	ArrayXXd distdiffmat = ArrayXXd::Zero(resn, resn), ref_distdiffmat = ArrayXXd::Zero(resn, resn);
	MatrixXd dist_holo = MatrixXd::Zero(resn, resn);
	// Coordinate gradient - 3N
	VectorXd gradient = VectorXd::Zero(veclen), off_pocket_gradient = VectorXd::Zero(veclen);
	// Equilibrium Coordinate - 3N
	VectorXd prev_equilibrium = VectorXd::Zero(veclen), ref_equilibrium_coord = VectorXd::Zero(veclen);

	flag = false;
	
	// Iteration
	for (size_t k = 0; k < paras.niteration; ++k)
	{
		dist_holo = gen_distmat(Dist, equilibrium_coord);
		distdiffmat = dist_holo - dist_apo;

		potential = (distdiffmat.pow(2) * kmat).sum() / 4;
		handle_message(MSG_EMPTY, boost::format("potential: %1$.4f") % potential);

		if (k > 0 && abs(potential - prev_potential) < paras.convergence)
		{
			flag = true;
			// Print gradient (vector)
			break;
		}

		prev_potential = potential;
		gradient = coords_derivate(equilibrium_coord, dist_apo, dist_holo, kmat);

		off_pocket_gradient = gradient;
		modify_pocket_coord(off_pocket_gradient, VectorXd::Zero(pocketlen), pocket_members);

		err = sqrt(off_pocket_gradient.dot(off_pocket_gradient));
		handle_message(MSG_RESULT, boost::format("Error: $1.6f") % err);

		ref_equilibrium_coord = equilibrium_coord - step * gradient / err;
		modify_pocket_coord(ref_equilibrium_coord, fix_coord, pocket_members);

		ref_distdiffmat = gen_distmat(Dist, ref_equilibrium_coord) - dist_apo;
		ref_potential = (ref_distdiffmat.pow(2) * kmat).sum() / 4;

		if (ref_potential < potential)
		{
			step *= 1.2;
			equilibrium_coord = ref_equilibrium_coord;
		}
		else if (ref_potential > potential)
		{
			step /= 1.5;
			prev_potential = 0;
		}
	}

	if (!flag)
	{
		handle_message(MSG_WARNING, boost::format("Do not converge in %1% steps.") % paras.niteration);
		return;
	}

	double rmsd_full = calc_rmsd(equilibrium_coord, fitting(initial_coord, equilibrium_coord));
	handle_message(MSG_RESULT, boost::format("RMSD for fitting full coordinates before and after minimization: $1.4f ") % rmsd_full);

	equilibrium_coord = fitting(initial_coord, equilibrium_coord);

	VectorXd pocket = VectorXd::Zero(pocketlen);
	VectorXd apo_pocket = VectorXd::Zero(pocketlen);
	grep_pocket_coord(pocket, equilibrium_coord, pocket_members);
	grep_pocket_coord(apo_pocket, initial_coord, pocket_members);

	double rmsd_pocket = calc_rmsd(pocket, fitting(apo_pocket, pocket));
	handle_message(MSG_RESULT, boost::format("RMSD for fitting pocket coordinates before and after minimization: $1.4f ") % rmsd_pocket);

	gradient = coords_derivate(equilibrium_coord, dist_apo, dist_holo, kmat);

	copy_pocket_coord(pocket_force, gradient, pocket_members);
}

void ProAnalysis::minimization(bool & flag, VectorXd & equilibrium_coord, VectorXd & pocket_force, BGDpara paras)
{
	double step = paras.learning_rate;

	// Set variables that would be often used
	VectorXd initial_coord = ProE.get_procoord();
	size_t veclen = initial_coord.size(); // 3N
	size_t resn = veclen / 3; // N
	ArrayXXd kmat = ProE.get_kmat();
	MatrixXd dist_apo = gen_distmat(Dist, initial_coord);

	equilibrium_coord = initial_coord;

	// Initialize variables
	double prev_potential = 0.0, potential = 0.0, ref_potential = 0.0, err = 0.0;
	// Pocket Force - 3N
	pocket_force = VectorXd::Zero(veclen);
	// Distance Matrix: Dist - N x N
	ArrayXXd distdiffmat = ArrayXXd::Zero(resn, resn), ref_distdiffmat = ArrayXXd::Zero(resn, resn);
	MatrixXd dist_holo = MatrixXd::Zero(resn, resn);
	// Coordinate gradient - 3N
	VectorXd gradient = VectorXd::Zero(veclen);
	// Equilibrium Coordinate - 3N
	VectorXd prev_equilibrium = VectorXd::Zero(veclen), ref_equilibrium_coord = VectorXd::Zero(veclen);

	flag = false;

	for (size_t k = 0; k < paras.niteration; ++k)
	{
		dist_holo = gen_distmat(Dist, equilibrium_coord);
		distdiffmat = dist_holo - dist_apo;

		potential = (distdiffmat.pow(2) * kmat).sum() / 4 - pocket_force.transpose() * (equilibrium_coord - initial_coord);
		handle_message(MSG_EMPTY, boost::format("potential: %1$.4f") % potential);

		if (k > 0 && abs(potential - prev_potential) < paras.convergence)
		{
			flag = true;
			// Print gradient (vector)
			break;
		}

		prev_potential = potential;

		gradient = coords_derivate(equilibrium_coord, dist_apo, dist_holo, kmat);		
		gradient -= pocket_force;

		err = sqrt(gradient.dot(gradient));

		ref_equilibrium_coord = equilibrium_coord - step * gradient / err;
		// Print step (double)

		ref_distdiffmat = gen_distmat(Dist, ref_equilibrium_coord) - dist_apo;
		double new_potential = (ref_distdiffmat.pow(2) * kmat).sum() / 4 - pocket_force.transpose() * (ref_equilibrium_coord - initial_coord);

		if (new_potential < potential)
		{
			step *= 1.2;
			equilibrium_coord = ref_equilibrium_coord;
		}
		else if (new_potential > potential)
		{
			step /= 1.5;
			prev_potential = 0;
		}
	}

	if (!flag)
	{
		handle_message(MSG_WARNING, boost::format("Do not converge in %1% steps.") % paras.niteration);
		return;
	}

	double rmsd_full = calc_rmsd(equilibrium_coord, fitting(initial_coord, equilibrium_coord));
	handle_message(MSG_RESULT, boost::format("RMSD for fitting full coordinates before and after minimization: $1.4f ") % rmsd_full);

	equilibrium_coord = fitting(initial_coord, equilibrium_coord);
}

void ProAnalysis::optimize_pocket_structure(Vector3d &vcenter, Vector3d &degrees, VectorXd coord, MatrixXd distmat_0, MatrixXd distmat, ArrayXXd kmat)
{
	assert(distmat_0.cols() == distmat_0.rows());
	assert(distmat.cols() == distmat.rows());
	assert(kmat.cols() == kmat.rows());
	assert(distmat_0.rows() == distmat.rows() == kmat.rows());
	assert(coord.size() == 3 * distmat_0.rows());

	size_t lenN = distmat_0.rows();

	VectorXd dRvec = coords_derivate(coord, distmat_0, distmat, kmat);
	VectorXd dDeg = degrees_derivate(coord, dRvec, vcenter, degrees);

	
}