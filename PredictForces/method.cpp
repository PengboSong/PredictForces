#include "method.h"

VectorXd rotate(VectorXd coord, Vector3d axis, double angle)
{
	if (axis.dot(axis) != 1.0)
		axis /= sqrt(axis.dot(axis));

	Matrix3d rotmat = Matrix3d::Zero();
	rotmat(1, 0) = axis(2);
	rotmat(0, 2) = axis(1);
	rotmat(2, 1) = axis(0);
	rotmat(0, 1) = -rotmat(1, 0);
	rotmat(2, 0) = -rotmat(0, 2);
	rotmat(1, 2) = -rotmat(2, 1);

	rotmat *= sin(angle);
	rotmat += (1 - cos(angle)) * (axis * axis.transpose());
	rotmat += cos(angle) * Matrix3d::Identity();

	size_t lenN = coord.size() / 3;
	Map<Matrix3Xd> xyz(coord.data(), 3, lenN);
	Vector3d center = xyz.rowwise().mean();
	xyz.colwise() -= center;

	Matrix3Xd rot_xyz = rotmat * xyz;
	rot_xyz.colwise() += center;
	Map<VectorXd> rot_coord(rot_xyz.data(), rot_xyz.size());

	return rot_coord;
}

VectorXd fitting(VectorXd coord1, VectorXd coord2)
{
	if (coord1.size() == coord2.size())
	{
		size_t vlen = coord1.size() / 3;
		Map<Matrix3Xd> xyz1(coord1.data(), 3, vlen);
		Map<Matrix3Xd> xyz2(coord2.data(), 3, vlen);
		Vector3d c1 = xyz1.rowwise().mean();
		Vector3d c2 = xyz2.rowwise().mean();
		xyz1.colwise() -= c1;
		xyz2.colwise() -= c2;

		Matrix3d r_mat = Matrix3d::Zero();
		for (size_t i = 0; i < 3; ++i)
			for (size_t j = 0; j < 3; ++j)
				for (size_t k = 0; k < vlen; ++k)
					r_mat(i, j) += xyz2(i, k) * xyz1(j, k);

		Matrix3d sym_mat = r_mat.transpose() * r_mat;

		SelfAdjointEigenSolver<Matrix3d> eigensolver(sym_mat);
		Vector3d eigenvalues = eigensolver.eigenvalues();
		Matrix3d eigenvectors = eigensolver.eigenvectors();

		Matrix3d b_mat = Matrix3d::Zero();
		for (size_t i = 0; i < 3; ++i)
			b_mat.row(i) = r_mat * eigenvectors.col(i) / sqrt(eigenvalues[i]);

		Matrix3d a_mat = eigenvectors.transpose();

		Matrix3d u_mat = Matrix3d::Zero();
		for (size_t i = 0; i < 3; ++i)
			for (size_t j = 0; j < 3; ++j)
				for (size_t k = 0; k < 3; ++k)
					u_mat(i, j) += b_mat(k, j) * a_mat(k, i);

		Matrix3Xd fit_xyz2 = u_mat * xyz2;
		fit_xyz2.colwise() += c1;
		Map<VectorXd> fit_coord2(fit_xyz2.data(), fit_xyz2.size());
		return fit_coord2;
	}
	else
		return VectorXd();
}

VectorXd calc_displacement(VectorXd coord1, VectorXd coord2)
{
	if (coord1.size() == coord2.size())
		return coord2 - coord1;
	else
		return VectorXd();
}

double calc_rmsd(VectorXd coord1, VectorXd coord2)
{
	if (coord1.size() == coord2.size())
	{
		VectorXd diff = coord1 - coord2;
		size_t vlen = diff.size() / 3;

		Map<Array3Xd> diffxyz(diff.data(), 3, vlen);

		VectorXd sqdist = diffxyz.pow(2).colwise().sum();
		return sqrt(sqdist.sum() / vlen);
	}
	else
		return 0.0;
}

double calc_rmsd(VectorXd diff)
{
	size_t vlen = diff.size() / 3;

	Map<Array3Xd> diffxyz(diff.data(), 3, vlen);

	VectorXd sqdist = diffxyz.pow(2).colwise().sum();
	return sqrt(sqdist.sum() / vlen);
}

double calc_mindist(VectorXd coord1, VectorXd coord2)
{
	Map<Matrix3Xd> xyz1(coord1.data(), 3, coord1.size() / 3);
	Map<Matrix3Xd> xyz2(coord2.data(), 3, coord2.size() / 3);

	VectorXd dist(xyz1.cols());
	Vector3d onexyz = Vector3d::Zero();
	Array3Xd diffxyz(3, xyz2.cols());
	for (size_t i = 0; i < size_t(xyz1.cols()); ++i)
	{
		onexyz = xyz1.col(i);
		diffxyz = xyz2.colwise() - onexyz;
		dist(i) = diffxyz.pow(2).colwise().sum().minCoeff();
	}
	return dist.minCoeff();
}

double calc_norm(VectorXd vector)
{
	ArrayXd v = vector;
	return sqrt(v.pow(2).sum());
}

double calc_average(VectorXd force)
{
	Map<Array3Xd> forcexyz(force.data(), 3, force.size() / 3);
	VectorXd sqdist = forcexyz.pow(2).colwise().sum().sqrt();
	return sqdist.sum() / sqdist.size();
}

MatrixXd gen_distmat(DistMatType t, VectorXd coord, double cutoff)
{
	size_t resn = coord.size() / 3;
	MatrixXd distmat;
	Vector3d coordi;
	Array3d coordij;
	double dist = 0.0;
	if (t == Dist)
	{
		distmat = MatrixXd::Zero(resn, resn);
		for (size_t i = 0; i < resn; ++i)
		{
			coordi = coord.block(3 * i, 1, 3, 1);
			for (size_t j = i + 1; j < resn; ++j)
			{
				coordij = coordi - coord.block(3 * j, 1, 3, 1);
				dist = sqrt(coordij.pow(2).sum());
				distmat(j, i) = distmat(i, j) = dist;
			}
		}
	}
	else if (t == XYZdiff || t == XYZdiffLtCutoff)
	{
		distmat = MatrixXd::Zero(resn, 3 * resn);
		if (t == XYZdiff)
		{
			for (size_t i = 0; i < resn; ++i)
			{
				coordi = coord.block(3 * i, 1, 3, 1);
				for (size_t j = i + 1; j < resn; ++j)
				{
					coordij = coordi - coord.block(3 * j, 1, 3, 1);
					distmat(i, 3 * j) = coordij(0);
					distmat(i, 3 * j + 1) = coordij(1);
					distmat(i, 3 * j + 2) = coordij(2);
					distmat(j, 3 * i) = -coordij(0);
					distmat(j, 3 * i + 1) = -coordij(1);
					distmat(j, 3 * i + 2) = -coordij(2);
				}
			}
		}
		else
		{
			for (size_t i = 0; i < resn; ++i)
			{
				coordi = coord.block(3 * i, 1, 3, 1);
				for (size_t j = i + 1; j < resn; ++j)
				{
					coordij = coordi - coord.block(3 * j, 1, 3, 1);
					dist = sqrt(coordij.pow(2).sum());
					if (dist < cutoff)
					{
						distmat(i, 3 * j) = coordij(0);
						distmat(i, 3 * j + 1) = coordij(1);
						distmat(i, 3 * j + 2) = coordij(2);
						distmat(j, 3 * i) = -coordij(0);
						distmat(j, 3 * i + 1) = -coordij(1);
						distmat(j, 3 * i + 2) = -coordij(2);
					}
				}
			}
		}
	}
	return distmat;
}

MatrixXd gen_differ(MatrixXd X, MatrixXd Y)
{
	if (X.cols() == Y.cols() && X.rows() == Y.rows())
	{
		MatrixXd diffmat = MatrixXd::Zero(X.rows(), X.cols());
		for (size_t i = 0; i < size_t(diffmat.rows()); i++)
			for (size_t j = i; j < size_t(diffmat.cols()); j++)
			{
				if (Y(i, j) != 0)
					diffmat(i, j) = X(i, j) - Y(i, j);
			}
		return diffmat;
	}
	else
	{
		handle_message(
			MSG_WARNING,
			"Can only generate a diff matrix for two matrixes with the same shape."
		);
		return MatrixXd();
	}
}

std::list<size_t> gen_pocket(double cutoff, VectorXd dist2ligand)
{
	std::list<size_t> pocket;
	if (cutoff > dist2ligand.minCoeff())
	{
		for (size_t i = 0; i < size_t(dist2ligand.size()); ++i)
			if (dist2ligand(i) < cutoff)
				pocket.push_back(i);
	}
	else
		handle_message(
			MSG_ERROR,
			boost::format("Given cutoff is too short. Minimum possible cutoff is %1$.2f.") % dist2ligand.minCoeff()
		);
	return pocket;
}

Matrix3d euler_rotation_matrix(Vector3d degrees)
{
	Matrix3d rmat = Matrix3d::Zero();

	double varphi = degrees[0], phi = degrees[1], theta = degrees[2];

	double cos_varphi = cos(varphi), sin_varphi = sin(varphi);
	double cos_phi = cos(phi), sin_phi = sin(phi);
	double cos_theta = cos(theta), sin_theta = sin(theta);

	rmat(0, 0) = cos_theta * cos_phi;
	rmat(1, 0) = cos_theta * sin_phi;
	rmat(2, 0) = -sin_theta;
	rmat(1, 0) = sin_varphi * sin_theta * cos_phi - cos_varphi * sin_phi;
	rmat(1, 1) = sin_varphi * sin_theta * sin_phi + cos_varphi * cos_phi;
	rmat(1, 2) = sin_varphi * cos_theta;
	rmat(2, 0) = cos_varphi * sin_theta * cos_phi + sin_varphi * sin_phi;
	rmat(2, 1) = cos_varphi * sin_theta * sin_phi - sin_varphi * cos_phi;
	rmat(2, 2) = cos_varphi * cos_theta;

	return rmat;
}

VectorXd coords_derivate(VectorXd coord, MatrixXd distmat_0, MatrixXd distmat, ArrayXXd kmat)
{
	assert(distmat_0.cols() == distmat_0.rows());
	assert(distmat.cols() == distmat.rows());
	assert(kmat.cols() == kmat.rows());
	assert(distmat_0.rows() == distmat.rows() == kmat.rows());
	assert(coord.size() == 3 * distmat_0.rows());

	size_t lenN = distmat_0.rows(); // N

	VectorXd dR = VectorXd::Zero(3 * lenN);

	// Distance Matrix : XYZdiff - N x 3N
	MatrixXd xyzdiff = gen_distmat(XYZdiff, coord);
	// Matrix { (dist - dist0) * k / dist0 } - N x N
	MatrixXd coeff = (distmat - distmat_0).array() * kmat / distmat.array();

	double coeffij = 0.0;
	for (size_t i = 0; i < lenN; ++i)
	{
		for (size_t j = i + 1; j < lenN; ++j)
		{
			coeffij = coeff(i, j);
			dR(3 * i) += xyzdiff(i, 3 * j) * coeffij;
			dR(3 * i + 1) += xyzdiff(i, 3 * j + 1) * coeffij;
			dR(3 * i + 2) += xyzdiff(i, 3 * j + 2) * coeffij;
		}
	}
	return dR;
}

Vector3d euler_degrees_derivate(Vector3d r, Vector3d dvdr, Vector3d degrees)
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

VectorXd degrees_derivate(VectorXd coord, VectorXd dR, Vector3d vcenter, Vector3d degrees)
{
	size_t lenN = coord.size() / 3;
	VectorXd dDeg = VectorXd::Zero(3 * lenN);
	for (size_t i = 0; i < lenN; ++i)
		dDeg.block(3 * i, 1, 3, 1) = euler_degrees_derivate(coord.block(3 * i, 1, 3, 1) - vcenter, dR.col(i), degrees);
	return dDeg;
}


void normal_equation(VectorXd &coeff, MatrixXd X, VectorXd Y)
{
	coeff = (X.transpose() * X).inverse() * X.transpose() * Y;
}

void BGD(VectorXd &coeff, MatrixXd X, VectorXd Y, BGDpara paras, bool checkinf)
{
	size_t nfeature = coeff.size();
	size_t nsample = Y.size();

	double step = paras.learning_rate;

	double cost = 0.0, prev_cost = 0.0;

	VectorXd gradient = VectorXd::Zero(nfeature);
	VectorXd new_gradient = VectorXd::Zero(nfeature);
	double gradient_product = 0.0, new_gradient_product = 0.0;

	bool converge_flag = false;

	for (size_t k = 0; k < paras.niteration; ++k)
	{
		handle_message(MSG_EMPTY, "*---*---*---*---*---*");
		handle_message(
			MSG_EMPTY,
			boost::format("Step: %1d") % (k + 1)
		);

		cost = (X * coeff - Y).dot(X * coeff - Y) / 2 / nsample;

		handle_message(
			MSG_EMPTY,
			boost::format("Cost: %1$.4f") % cost
		);

		if (k > 0 && abs(cost - prev_cost) < paras.convergence)
		{
			converge_flag = true;
			handle_message(
				MSG_EMPTY,
				boost::format("Gradient: %1$.6f") % gradient
			);
			break;
		}

		handle_message(MSG_EMPTY, "*---*---*---*---*---*");

		prev_cost = cost;
		gradient = ((X * coeff - Y).transpose() * X).transpose() / nsample;
		gradient_product = gradient.dot(gradient);
		coeff -= step / nsample * ((X * coeff - Y).transpose() * X).transpose();
		new_gradient = ((X * coeff - Y).transpose() * X).transpose() / nsample;
		new_gradient_product = new_gradient.dot(new_gradient);

		if (new_gradient_product >= gradient_product)
			step /= 2;
		else if (new_gradient_product < gradient_product)
			step *= 1.5;

		if (checkinf)
		{
			if (isinf<double>(cost))
			{
				handle_message(
					MSG_WARNING,
					boost::format("Reach infinity in %1% steps.") % (k + 1)
				);
				break;
			}
		}
	}

	if (!converge_flag)
	{
		handle_message(
			MSG_WARNING,
			boost::format("BGD does not converge in %1% steps.") % paras.niteration
		);
	}
}
