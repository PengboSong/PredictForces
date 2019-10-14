#include "ProAnalysis.h"

void ProAnalysis::minimization(bool & flag, PocketList pocket_members, VectorXd fix_coord, VectorXd & equilibrium_coord, VectorXd & pocket_force, BGDpara paras)
{
	double step = paras.learning_rate;

	VectorXd initial_coord = ProE.get_procoord();
	equilibrium_coord = initial_coord;
	size_t i = 0;
	for (PocketList::iterator it = pocket_members.begin(); it != pocket_members.end(); ++it)
	{
		equilibrium_coord(*it * 3) = fix_coord(i * 3);
		equilibrium_coord(*it * 3 + 1) = fix_coord(i * 3 + 1);
		equilibrium_coord(*it * 3 + 2) = fix_coord(i * 3 + 2);
		++i;
	}
	size_t resn = initial_coord.size();
	pocket_force = VectorXd::Zero(resn);
	double prev_potential = 0.0, potential = 0.0;
	double err = 0;
	bool converge_flag = false;
	ArrayXXd distdiffmat = ArrayXXd::Zero(resn, resn);
	MatrixXd apo = gen_small_distmat(initial_coord);
	ArrayXXd kmat = ProE.get_kmat();
	VectorXd gradient = VectorXd::Zero(resn);
	VectorXd prev_equilibrium = VectorXd::Zero(resn);

	for (size_t k = 0; k < paras.niteration; ++k)
	{
		MatrixXd holo = gen_all_distmat(equilibrium_coord);
		MatrixXd small_holo = gen_small_distmat(equilibrium_coord);
		distdiffmat = small_holo - apo;

		potential = (distdiffmat.pow(2) * kmat).sum() / 4;
		handle_message(MSG_EMPTY, boost::format("potential: %1$.4f") % potential);

		if (k > 0 && abs(potential - prev_potential) < paras.convergence)
		{
			converge_flag = true;
			std::cout << "Gradient:" << gradient << std::endl;
			break;
		}

		prev_potential = potential;
		gradient = VectorXd::Zero(resn);
		MatrixXd C = distdiffmat * kmat / small_holo.array();

		for (size_t i = 0; (i < resn / 3); ++i)
			for (size_t j = 0; (j < resn / 3); ++j)
				if (i != j)
				{
					gradient(3 * i) += holo(3 * i, 3 * j) * C(i, j);
					gradient(3 * i + 1) += holo(3 * i + 1, 3 * j + 1) * C(i, j);
					gradient(3 * i + 2) += holo(3 * i + 2, 3 * j + 2) * C(i, j);
				}
		VectorXd E = gradient;
		for (PocketList::iterator it = pocket_members.begin(); it != pocket_members.end(); ++it)
		{
			E(*it * 3) = 0;
			E(*it * 3 + 1) = 0;
			E(*it * 3 + 2) = 0;
		}
		err = sqrt(E.dot(E));
		std::cout << err << std::endl;

		VectorXd new_equilibrium_coord = equilibrium_coord - step * gradient / err;

		size_t j = 0;
		for (PocketList::iterator it = pocket_members.begin(); it != pocket_members.end(); ++it)
		{
			new_equilibrium_coord(*it * 3) = fix_coord(j * 3);
			new_equilibrium_coord(*it * 3 + 1) = fix_coord(j * 3 + 1);
			new_equilibrium_coord(*it * 3 + 2) = fix_coord(j * 3 + 2);
			++j;
		}

		MatrixXd new_small_holo = gen_small_distmat(new_equilibrium_coord);
		ArrayXXd new_distdiffmat = new_small_holo - apo;
		double new_potential = (new_distdiffmat.pow(2) * kmat).sum() / 4;

		if (new_potential < potential)
		{
			step *= 1.2;
			equilibrium_coord = new_equilibrium_coord;
		}
		else if (new_potential > potential)
		{
			step /= 1.5;
			prev_potential = 0;
		}
	}
	if (!converge_flag)
		handle_message(MSG_WARNING, boost::format("Do not converge in %1% steps.") % paras.niteration);

	double gall = calc_rmsd(equilibrium_coord, fitting(initial_coord, equilibrium_coord));
	std::cout << gall << std::endl;

	equilibrium_coord = fitting(initial_coord, equilibrium_coord);
	MatrixXd holo = gen_all_distmat(equilibrium_coord);
	MatrixXd small_holo = gen_small_distmat(equilibrium_coord);
	distdiffmat = small_holo - apo;
	gradient = VectorXd::Zero(resn);
	MatrixXd C = distdiffmat * kmat / small_holo.array();

	VectorXd poc = VectorXd::Zero(pocket_members.size() * 3);
	VectorXd apo_poc = VectorXd::Zero(pocket_members.size() * 3);
	size_t z = 0;
	for (PocketList::iterator it = pocket_members.begin(); it != pocket_members.end(); ++it)
	{
		poc(z * 3) = equilibrium_coord(*it * 3);
		poc(z * 3 + 1) = equilibrium_coord(*it * 3 + 1);
		poc(z * 3 + 2) = equilibrium_coord(*it * 3 + 2);

		apo_poc(z * 3) = initial_coord(*it * 3);
		apo_poc(z * 3 + 1) = initial_coord(*it * 3 + 1);
		apo_poc(z * 3 + 2) = initial_coord(*it * 3 + 2);
		++z;
	}
	double iall = calc_rmsd(poc, fitting(apo_poc, poc));
	std::cout << iall << std::endl;

	for (size_t i = 0; (i < resn / 3); ++i)
		for (size_t j = 0; (j < resn / 3); ++j)
			if (i != j)
			{
				gradient(3 * i) += holo(3 * i, 3 * j) * C(i, j);
				gradient(3 * i + 1) += holo(3 * i + 1, 3 * j + 1) * C(i, j);
				gradient(3 * i + 2) += holo(3 * i + 2, 3 * j + 2) * C(i, j);
			}

	for (PocketList::iterator it = pocket_members.begin(); it != pocket_members.end(); ++it)
	{
		pocket_force(*it * 3) = gradient(*it * 3);
		pocket_force(*it * 3 + 1) = gradient(*it * 3 + 1);
		pocket_force(*it * 3 + 2) = gradient(*it * 3 + 2);
	}

	flag = true;
}

void ProAnalysis::minimization(bool & flag, VectorXd & equilibrium_coord, VectorXd & pocket_force, BGDpara paras)
{
	double step = paras.learning_rate;

	VectorXd initial_coord = ProE.get_procoord();
	equilibrium_coord = initial_coord;

	size_t resn = initial_coord.size();
	VectorXd gradient = VectorXd::Zero(resn);
	double prev_potential = 0, potential = 0.0;
	double err = 0.0;
	bool converge_flag = false;
	MatrixXd apo = gen_small_distmat(initial_coord);
	ArrayXXd distdiffmat = ArrayXXd::Zero(resn, resn);
	ArrayXXd kmat = ProE.get_kmat();
	VectorXd prev_equilibrium = VectorXd::Zero(resn);

	for (size_t k = 0; k < paras.niteration; ++k)
	{
		MatrixXd holo = gen_all_distmat(equilibrium_coord);
		MatrixXd small_holo = gen_small_distmat(equilibrium_coord);
		distdiffmat = small_holo - apo;
		VectorXd displacement = equilibrium_coord - initial_coord;

		potential = (distdiffmat.pow(2) * kmat).sum() / 4 - pocket_force.transpose() * displacement;
		handle_message(MSG_EMPTY, boost::format("potential: %1$.4f") % potential);

		if (k > 0 && abs(potential - prev_potential) < paras.convergence)
		{
			converge_flag = true;
			std::cout << "Gradient:" << gradient << std::endl;
			break;
		}

		prev_potential = potential;
		gradient = VectorXd::Zero(resn);
		MatrixXd C = distdiffmat * kmat / small_holo.array();

		for (size_t i = 0; i < resn / 3; ++i)
			for (size_t j = 0; (j < resn / 3); ++j)
				if (i != j)
				{
					gradient(3 * i) += holo(3 * i, 3 * j) * C(i, j);
					gradient(3 * i + 1) += holo(3 * i + 1, 3 * j + 1) * C(i, j);
					gradient(3 * i + 2) += holo(3 * i + 2, 3 * j + 2) * C(i, j);
				}

		gradient -= pocket_force;
		err = sqrt(gradient.dot(gradient));

		VectorXd new_equilibrium_coord = equilibrium_coord - step * gradient / err;
		std::cout << step << std::endl;

		MatrixXd new_small_holo = gen_small_distmat(new_equilibrium_coord);
		ArrayXXd new_distdiffmat = new_small_holo - apo;
		VectorXd new_displacement = new_equilibrium_coord - initial_coord;
		double new_potential = (new_distdiffmat.pow(2) * kmat).sum() / 4 - pocket_force.transpose() * new_displacement;

		if (new_potential < potential)
		{
			step *= 1.2;
			equilibrium_coord = new_equilibrium_coord;
		}
		else if (new_potential > potential)
		{
			step /= 1.5;
			prev_potential = 0;
		}
	}
	if (!converge_flag)
		handle_message(MSG_WARNING, boost::format("Do not converge in %1% steps.") % paras.niteration);

	double gall = calc_rmsd(equilibrium_coord, fitting(initial_coord, equilibrium_coord));
	std::cout << gall << std::endl;
	equilibrium_coord = fitting(initial_coord, equilibrium_coord);

	flag = true;
}

void ProAnalysis::optimize_pocket_structure(Vector3d &vcenter, Vector3d &degrees, Matrix3Xd coord, MatrixXd distmat_0, MatrixXd distmat, MatrixXi contactmap, double k)
{
	Matrix3Xd dR = coords_derivate(coord, distmat_0, distmat, contactmap, k);

	Matrix3Xd dDeg = Matrix3Xd::Zero(3, dR.cols());
	for (size_t i = 0; i < size_t(dR.cols()); ++i)
	{
		dDeg.col(i) = euler_degrees_derivative(
			coord.col(i) - vcenter,
			dR.col(i),
			degrees
		);
	}
}