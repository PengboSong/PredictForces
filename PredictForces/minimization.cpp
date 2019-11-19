#include "Minimization.h"

Minimization::Minimization()
{
}

Minimization::Minimization(BGDpara paras, size_t N, VectorXd apo_procoord, ArrayXXd k_matrix)
{
	parameters = paras;
	noden = N;
	coordn = 3 * N;

	// Equilibrium Coordinate - 3N
	coord = VectorXd::Zero(coordn);
	coord_apo = apo_procoord;
	assert(coord_apo.size() == coordn);
	coord_ref = VectorXd::Zero(coordn);
	// Pocket Force - 3N
	pocketforce = VectorXd::Zero(coordn);
	// Distance Matrix: Dist - N x N
	distmat = MatrixXd::Zero(noden, noden);
	distmat_apo = gen_distmat(Dist, apo_procoord);
	assert(distmat_apo.rows() == noden && distmat_apo.cols() == noden);
	distmat_ref = MatrixXd::Zero(noden, noden);
	// Distance Matrix : XYZdiff - N x 3N
	xyzdiffmat = MatrixXd::Zero(noden, coordn);
	// Matrix { (dist - dist0) * k / dist0 } - N x N
	coeff = MatrixXd::Zero(noden, noden);

	kmat = k_matrix;
	assert(kmat.rows() == noden && kmat.cols() == noden);

	// Coordinate gradient - 3N
	dR = VectorXd::Zero(coordn);
}

Minimization::~Minimization()
{
}

void Minimization::update_coords_derivate()
{
	xyzdiffmat = gen_distmat(XYZdiff, coord);
	distmat = gen_distmat(Dist, coord);
	coeff = (distmat - distmat_apo).array() * kmat / distmat.array();

	dR = VectorXd::Zero(coordn);

	double coeffij = 0.0;
	for (size_t i = 0; i < noden; ++i)
	{
		for (size_t j = 0; j < noden; ++j)
		{
			if (i != j)
			{
				coeffij = coeff(i, j);
				dR(3 * i) += xyzdiffmat(i, 3 * j) * coeffij;
				dR(3 * i + 1) += xyzdiffmat(i, 3 * j + 1) * coeffij;
				dR(3 * i + 2) += xyzdiffmat(i, 3 * j + 2) * coeffij;
			}
		}
	}
}

Vector3d Minimization::euler_degrees_derivate(Vector3d r, Vector3d dvdr, Vector3d degrees)
{
	double rx = r[0], ry = r[1], rz = r[2];
	double dv_drx = dvdr[0], dv_dry = dvdr[1], dv_drz = dvdr[2];
	double varphi = degrees[0], phi = degrees[1], theta = degrees[2];

	double cos_varphi = cos(varphi), sin_varphi = sin(varphi);
	double cos_phi = cos(phi), sin_phi = sin(phi);
	double cos_theta = cos(theta), sin_theta = sin(theta);

	double varphi_cymsz = cos_varphi * ry - sin_varphi * rz;
	double varphi_syacz = sin_varphi * ry + cos_varphi * rz;

	double drx_dvarphi = sin_theta * cos_phi * varphi_cymsz + sin_phi * varphi_syacz;
	double dry_dvarphi = sin_theta * sin_phi * varphi_cymsz - cos_phi * varphi_syacz;
	double drz_dvarphi = cos_theta * varphi_cymsz;

	double drx_dphi = -cos_theta * sin_phi * rx - sin_theta * sin_phi * varphi_syacz - cos_phi * varphi_cymsz;
	double dry_dphi = cos_theta * cos_phi * rx + sin_theta * cos_phi * varphi_syacz - sin_phi * varphi_cymsz;
	double drz_dphi = 0.0;

	double drx_dtheta = -cos_phi * sin_theta * rx + cos_phi * cos_theta * varphi_syacz;
	double dry_dtheta = -sin_phi * sin_theta * rx + sin_phi * cos_theta * varphi_syacz;
	double drz_dtheta = -cos_theta * rx - sin_theta * varphi_syacz;

	double dv_dvarphi = dv_drx * drx_dvarphi + dv_dry * dry_dvarphi + dv_drz * drz_dvarphi;
	double dv_dphi = dv_drx * drx_dphi + dv_dry * dry_dphi + dv_drz * drz_dphi;
	double dv_dtheta = dv_drx * drx_dtheta + dv_dry * dry_dtheta + dv_drz * drz_dtheta;

	Vector3d dDegv;
	dDegv << dv_dvarphi, dv_dphi, dv_dtheta;
	return dDegv;
}

VectorXd Minimization::degrees_derivate(VectorXd domain_coord, Vector3d vcenter, Vector3d degrees)
{
	VectorXd dDeg = VectorXd::Zero(coordn);
	for (size_t i = 0; i < noden; ++i)
		dDeg.block(3 * i, 0, 3, 1) = euler_degrees_derivate(domain_coord.block(3 * i, 0, 3, 1) - vcenter, dR.col(i), degrees);
	return dDeg;
}

void Minimization::print_xyzvec(VectorXd vec, std::string prefix)
{
	assert(vec.size() == coordn);
	Log.handle_message(MSG_EMPTY, "*---*---*---*---*---*");
	for (size_t id = 0; id < noden; ++id)
		Log.handle_message(
			MSG_EMPTY,
			boost::format("%1% %2$4f: [%3$9.4f, %4$9.4f, %5$9.4f]") % prefix % id % dR(3 * id) % dR(3 * id + 1) % dR(3 * id + 2)
		);
	Log.handle_message(MSG_EMPTY, "*---*---*---*---*---*");
}

void Minimization::reset()
{
	coord = VectorXd::Zero(coordn);
	coord_ref = VectorXd::Zero(coordn);
	pocketforce = VectorXd::Zero(coordn);
	distmat = MatrixXd::Zero(noden, noden);
	distmat_ref = MatrixXd::Zero(noden, noden);
	xyzdiffmat = MatrixXd::Zero(noden, coordn);
	dR = VectorXd::Zero(coordn);

	realstep = parameters.learning_rate;

	convergeStatus = false;
	completeFlag = false;
	downwardStep = false;
}

void Minimization::single_pocket(PocketList pocket_members, VectorXd fix_coord)
{
	reset();

	coord = coord_apo;
	modify_pocket_coord(coord, fix_coord, pocket_members);

	// Pocket length - L
	size_t pocketlen = pocket_members.size() * 3;
	// Coordinate gradient - 3N
	VectorXd dR_pocketoff = VectorXd::Zero(coordn);

	Log.handle_message(MSG_INFO, "Start potential minimization process.");
	Log.handle_message(MSG_EMPTY, "*---*---*---*---*---*");

	// Iteration
	for (size_t k = 0; k < parameters.niteration; ++k)
	{
		Log.handle_message(MSG_EMPTY, boost::format("Step: %1%") % k);

		distmat = gen_distmat(Dist, coord);

		potential = ((distmat - distmat_apo).array().pow(2) * kmat).sum() / 4;
		Log.handle_message(MSG_EMPTY, boost::format("Potential: %1$.4f") % potential);

		if (k > 0 && !downwardStep && abs(potential - prev_potential) < parameters.convergence)
		{
			convergeStatus = true;
			Log.handle_message(MSG_RESULT, boost::format("Complete potential minimization process within %1% steps.") % k);
			// Print gradient (vector)
			Log.handle_message(MSG_RESULT, "Gradient of the last step:");
			print_gradient();
			break;
		}

		prev_potential = potential;
		update_coords_derivate();

		dR_pocketoff = dR;
		modify_pocket_coord(dR_pocketoff, VectorXd::Zero(pocketlen), pocket_members);

		err = sqrt(dR_pocketoff.dot(dR_pocketoff));
		Log.handle_message(MSG_EMPTY, boost::format("Error: %1$.4f") % err);

		coord_ref = coord - realstep * dR / err;
		modify_pocket_coord(coord_ref, fix_coord, pocket_members);
		Log.handle_message(MSG_EMPTY, boost::format("Step Length: %1$.3f") % realstep);
		Log.handle_message(MSG_EMPTY, "*---*---*---*---*---*");

		distmat_ref = gen_distmat(Dist, coord_ref);
		ref_potential = ((distmat_ref - distmat_apo).array().pow(2) * kmat).sum() / 4;

		if (ref_potential < potential)
		{
			realstep *= parameters.upward_factor;
			coord = coord_ref;
			downwardStep = false;
		}
		else if (ref_potential > potential)
		{
			realstep /= parameters.downward_factor;
			downwardStep = true;
		}
	}

	if (!convergeStatus)
	{
		Log.handle_message(MSG_WARNING, boost::format("Do not converge in %1% steps.") % parameters.niteration);
		return;
	}

	double rmsd_full = calc_rmsd(coord, fitting(coord_apo, coord));
	Log.handle_message(MSG_RESULT, boost::format("RMSD for fitting full coordinates before and after minimization: %1$.4f A.") % rmsd_full);

	coord = fitting(coord_apo, coord);
	Log.handle_message(MSG_RESULT, "Equilibrium coordinates:");
	print_coord();

	VectorXd pocket = VectorXd::Zero(pocketlen);
	VectorXd pocket_apo = VectorXd::Zero(pocketlen);
	grep_pocket_coord(pocket, coord, pocket_members);
	grep_pocket_coord(pocket_apo, coord_apo, pocket_members);

	double rmsd_pocket = calc_rmsd(pocket, fitting(pocket_apo, pocket));
	Log.handle_message(MSG_RESULT, boost::format("RMSD for fitting pocket coordinates before and after minimization: %1$.4f A.") % rmsd_pocket);

	update_coords_derivate();
	copy_pocket_coord(pocketforce, dR, pocket_members);
	Log.handle_message(MSG_RESULT, "Pocket force:");
	print_pocketforce();

	completeFlag = true;
}

void Minimization::single_pocket(VectorXd pocket_force)
{
	reset();

	coord = coord_apo;

	// Pocket Force - 3N
	assert(pocket_force.size() == coordn);
	pocketforce = pocket_force;
	// Equilibrium Coordinate - 3N
	VectorXd coord_ref = VectorXd::Zero(coordn);

	Log.handle_message(MSG_INFO, "Start potential minimization process.");
	Log.handle_message(MSG_EMPTY, "*---*---*---*---*---*");

	for (size_t k = 0; k < parameters.niteration; ++k)
	{
		Log.handle_message(MSG_EMPTY, boost::format("Step: %1%") % k);

		distmat = gen_distmat(Dist, coord);

		potential = ((distmat - distmat_apo).array().pow(2) * kmat).sum() / 4 - pocketforce.dot(coord - coord_apo);
		Log.handle_message(MSG_EMPTY, boost::format("Potential: %1$.4f") % potential);

		if (k > 0 && !downwardStep && abs(potential - prev_potential) < parameters.convergence)
		{
			convergeStatus = true;
			Log.handle_message(MSG_RESULT, boost::format("Complete potential minimization process within %1% steps.") % k);
			// Print gradient (vector)
			Log.handle_message(MSG_RESULT, "Gradient of the last step:");
			print_gradient();
			break;
		}

		prev_potential = potential;

		update_coords_derivate();
		dR -= pocketforce;

		err = sqrt(dR.dot(dR));
		Log.handle_message(MSG_EMPTY, boost::format("Error: %1$.4f") % err);

		coord_ref = coord - realstep * dR / err;
		Log.handle_message(MSG_EMPTY, boost::format("Step Length: %1$.3f") % realstep);
		Log.handle_message(MSG_EMPTY, "*---*---*---*---*---*");

		distmat_ref = gen_distmat(Dist, coord_ref);
		ref_potential = ((distmat_ref - distmat_apo).array().pow(2) * kmat).sum() / 4 - pocketforce.dot(coord_ref - coord_apo);

		if (ref_potential < potential)
		{
			realstep *= parameters.upward_factor;
			coord = coord_ref;
			downwardStep = false;
		}
		else if (ref_potential > potential)
		{
			realstep /= parameters.downward_factor;
			downwardStep = true;
		}
	}

	if (!convergeStatus)
	{
		Log.handle_message(MSG_WARNING, boost::format("Do not converge in %1% steps.") % parameters.niteration);
		return;
	}

	double rmsd_full = calc_rmsd(coord, fitting(coord_apo, coord));
	Log.handle_message(MSG_RESULT, boost::format("RMSD for fitting full coordinates before and after minimization: %1$.4f A.") % rmsd_full);

	coord = fitting(coord_apo, coord);
	Log.handle_message(MSG_RESULT, "Equilibrium coordinates:");
	print_coord();

	Log.handle_message(MSG_RESULT, "Pocket force:");
	print_pocketforce();

	completeFlag = true;
}

void Minimization::single_pocket(PocketList pocket_members, VectorXd pocket_force, VectorXd fix_coord)
{
	reset();
	
	coord = coord_apo;
	modify_pocket_coord(coord, fix_coord, pocket_members);

	// Pocket length - L
	size_t pocketlen = pocket_members.size() * 3;
	// Pocket Force - 3N
	assert(pocket_force.size() == coordn);
	pocketforce = pocket_force;
	// Coordinate gradient - 3N
	VectorXd dR_pocketoff = VectorXd::Zero(coordn);

	Log.handle_message(MSG_INFO, "Start potential minimization process.");
	Log.handle_message(MSG_EMPTY, "*---*---*---*---*---*");

	// Iteration
	for (size_t k = 0; k < parameters.niteration; ++k)
	{
		Log.handle_message(MSG_EMPTY, boost::format("Step: %1%") % k);

		distmat = gen_distmat(Dist, coord);

		potential = ((distmat - distmat_apo).array().pow(2) * kmat).sum() / 4 - pocketforce.dot(coord - coord_apo);
		Log.handle_message(MSG_EMPTY, boost::format("Potential: %1$.4f") % potential);

		if (k > 0 && !downwardStep && abs(potential - prev_potential) < parameters.convergence)
		{
			convergeStatus = true;
			Log.handle_message(MSG_RESULT, boost::format("Complete potential minimization process within %1% steps.") % k);
			// Print gradient (vector)
			Log.handle_message(MSG_RESULT, "Gradient of the last step:");
			print_gradient();
			break;
		}

		prev_potential = potential;
		update_coords_derivate();
		dR -= pocketforce;

		dR_pocketoff = dR;
		modify_pocket_coord(dR_pocketoff, VectorXd::Zero(pocketlen), pocket_members);
		err = sqrt(dR_pocketoff.dot(dR_pocketoff));
		Log.handle_message(MSG_EMPTY, boost::format("Error: %1$.4f") % err);

		coord_ref = coord - realstep * dR / err;
		modify_pocket_coord(coord_ref, fix_coord, pocket_members);
		Log.handle_message(MSG_EMPTY, boost::format("Step Length: %1$.3f") % realstep);
		Log.handle_message(MSG_EMPTY, "*---*---*---*---*---*");

		distmat_ref = gen_distmat(Dist, coord_ref);
		ref_potential = ((distmat_ref - distmat_apo).array().pow(2) * kmat).sum() / 4 - pocketforce.dot(coord_ref - coord_apo);

		if (ref_potential < potential)
		{
			realstep *= parameters.upward_factor;
			coord = coord_ref;
			downwardStep = false;
		}
		else if (ref_potential > potential)
		{
			realstep /= parameters.downward_factor;
			downwardStep = true;
		}
	}

	if (!convergeStatus)
	{
		Log.handle_message(MSG_WARNING, boost::format("Do not converge in %1% steps.") % parameters.niteration);
		return;
	}

	double rmsd_full = calc_rmsd(coord, fitting(coord_apo, coord));
	Log.handle_message(MSG_RESULT, boost::format("RMSD for fitting full coordinates before and after minimization: %1$.4f A.") % rmsd_full);

	coord = fitting(coord_apo, coord);
	Log.handle_message(MSG_RESULT, "Equilibrium coordinates:");
	print_coord();

	VectorXd pocket = VectorXd::Zero(pocketlen);
	VectorXd pocket_apo = VectorXd::Zero(pocketlen);
	grep_pocket_coord(pocket, coord, pocket_members);
	grep_pocket_coord(pocket_apo, coord_apo, pocket_members);

	double rmsd_pocket = calc_rmsd(pocket, fitting(pocket_apo, pocket));
	Log.handle_message(MSG_RESULT, boost::format("RMSD for fitting pocket coordinates before and after minimization: %1$.4f A.") % rmsd_pocket);

	VectorXd pocket_force_change = VectorXd::Zero(coordn);
	copy_pocket_coord(pocket_force_change, dR, pocket_members);
	pocketforce += pocket_force_change;
	Log.handle_message(MSG_RESULT, "Pocket force:");
	print_pocketforce();

	completeFlag = true;
}

void Minimization::multiple_pocket(
	PocketList fixed_pocket_members,
	VectorXd fixed_pocket_coord,
	VectorXd fixed_pocket_force,
	PocketList rotate_pocket_members,
	VectorXd rotate_pocket_coord,
	VectorXd rotate_pocket_force,
	Vector3d & vcenter,
	Vector3d & degrees)
{
	reset();

	size_t pocketn_fixed = fixed_pocket_members.size();
	assert(fixed_pocket_coord.size() == 3 * pocketn_fixed);
	size_t pocketn_rot = rotate_pocket_members.size();
	assert(rotate_pocket_coord.size() == 3 * pocketn_rot);
	Map<Matrix3Xd> rot_pocket_domain(rotate_pocket_coord.data(), 3, pocketn_rot);

	coord = coord_apo;
	modify_pocket_coord(coord, fixed_pocket_coord, fixed_pocket_members);
	modify_pocket_coord(coord, rotate_pocket_coord, rotate_pocket_members);

	Vector3d vcenter_rot = coord_center(rotate_pocket_coord);
	vcenter = vcenter_rot;
	degrees = Vector3d::Zero();

	// Derivates
	VectorXd dRdomain = VectorXd::Zero(3 * pocketn_rot), dDeg = VectorXd::Zero(3 * pocketn_rot);
	Vector3d dcenter = Vector3d::Zero(), ddegrees = Vector3d::Zero();

	Log.handle_message(MSG_INFO, "Start potential minimization process.");
	Log.handle_message(MSG_EMPTY, "*---*---*---*---*---*");

	for (size_t k = 0; k < parameters.niteration; ++k)
	{
		Log.handle_message(MSG_EMPTY, boost::format("Step: %1%") % k);

		distmat = gen_distmat(Dist, coord);

		potential = ((distmat - distmat_apo).array().pow(2) * kmat).sum() / 4;
		Log.handle_message(MSG_EMPTY, boost::format("potential: %1$.4f") % potential);

		if (k > 0 && !downwardStep && abs(potential - prev_potential) < parameters.convergence)
		{
			convergeStatus = true;
			// Print gradient (vector)
			Log.handle_message(MSG_RESULT, "Gradient of the last step:");
			print_gradient();
			break;
		}

		prev_potential = potential;

		update_coords_derivate();
		modify_pocket_coord(dR, VectorXd::Zero(3 * pocketn_fixed), fixed_pocket_members);
		modify_pocket_coord(dR, VectorXd::Zero(3 * pocketn_rot), rotate_pocket_members);

		grep_pocket_coord(dRdomain, dR, rotate_pocket_members);
		dDeg = degrees_derivate(rot_pocket_domain, vcenter, degrees);
		dcenter = coord_vec2mat(dRdomain).rowwise().sum();
		ddegrees = coord_vec2mat(dDeg).rowwise().sum();

		err = sqrt(dR.dot(dR));
		Log.handle_message(MSG_EMPTY, boost::format("Error: %1$.4f") % err);

		rot_pocket_domain.colwise() -= vcenter;
		vcenter -= realstep * dcenter / err;
		degrees -= realstep * dDeg / err;
		rot_pocket_domain = euler_rotation_matrix(degrees) * rot_pocket_domain;
		rot_pocket_domain.colwise() += vcenter;

		coord_ref = coord - realstep * dR / err;
		modify_pocket_coord(coord_ref, coord_mat2vec(rot_pocket_domain), rotate_pocket_members);
		Log.handle_message(MSG_EMPTY, boost::format("Step Length: %1$.3f") % realstep);
		Log.handle_message(MSG_EMPTY, "*---*---*---*---*---*");

		distmat_ref = gen_distmat(Dist, coord_ref);
		ref_potential = ((distmat_ref - distmat_apo).array().pow(2) * kmat).sum() / 4;

		if (ref_potential < potential)
		{
			realstep *= 1.2;
			coord = coord_ref;
			downwardStep = false;
		}
		else if (ref_potential > potential)
		{
			realstep /= 1.5;
			downwardStep = true;
		}
	}

	if (!convergeStatus)
	{
		Log.handle_message(MSG_WARNING, boost::format("Do not converge in %1% steps.") % parameters.niteration);
		return;
	}

	double rmsd_full = calc_rmsd(coord, fitting(coord_apo, coord));
	Log.handle_message(MSG_RESULT, boost::format("RMSD for fitting full coordinates before and after minimization: %1$.4f ") % rmsd_full);

	coord = fitting(coord_apo, coord);

	completeFlag = true;
}

FreeEnergy Minimization::energy()
{
	FreeEnergy G;
	if (completeFlag)
	{
		distmat = gen_distmat(Dist, coord);
		G.pro = ((distmat - distmat_apo).array().pow(2) * kmat).sum() / 4;
		G.pocket = -pocketforce.dot(coord - coord_apo);
		G.total = G.pro + G.pocket;
	}
	else
	{
		Log.handle_message(MSG_WARNING, "No equilibrium coordinates from energy minimization task. Perform an energy minimization task before calling this method.");
		G.pro = G.pocket = G.total = 0.0;
	}
	return G;
}
