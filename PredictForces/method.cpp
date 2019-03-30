#include "method.h"

Eigen::VectorXd rotate(Eigen::VectorXd coord, Eigen::Vector3d axis, double angle)
{
	if (axis.dot(axis) != 1.0)
		axis /= sqrt(axis.dot(axis));

	Eigen::Matrix3d rotmat = Eigen::Matrix3d::Zero();
	rotmat(1, 0) = axis(2);
	rotmat(0, 2) = axis(1);
	rotmat(2, 1) = axis(0);
	rotmat(0, 1) = -rotmat(1, 0);
	rotmat(2, 0) = -rotmat(0, 2);
	rotmat(1, 2) = -rotmat(2, 1);

	rotmat *= sin(angle);
	rotmat += (1 - cos(angle)) * (axis * axis.transpose());
	rotmat += cos(angle) * Eigen::Matrix3d::Identity();
	
	Eigen::Map<Eigen::Matrix3Xd> xyz(coord.data(), 3, coord.size() / 3);
	Eigen::Vector3d center = xyz.rowwise().mean();
	xyz.colwise() -= center;

	Eigen::Matrix3Xd rot_xyz = rotmat * xyz;
	rot_xyz.colwise() += center;
	Eigen::Map<Eigen::VectorXd> rot_coord(rot_xyz.data(), rot_xyz.size());

	return rot_coord;
}

Eigen::VectorXd fitting(Eigen::VectorXd coord1, Eigen::VectorXd coord2)
{
	if (coord1.size() == coord2.size())
	{
		size_t vlen = coord1.size() / 3;
		Eigen::Map<Eigen::Matrix3Xd> xyz1(coord1.data(), 3, vlen);
		Eigen::Map<Eigen::Matrix3Xd> xyz2(coord2.data(), 3, vlen);
		Eigen::Vector3d c1 = xyz1.rowwise().mean();
		Eigen::Vector3d c2 = xyz2.rowwise().mean();
		xyz1.colwise() -= c1;
		xyz2.colwise() -= c2;

		Eigen::Matrix3d r_mat = Eigen::Matrix3d::Zero();
		for (size_t i = 0; i < 3; ++i)
			for (size_t j = 0; j < 3; ++j)
				for (size_t k = 0; k < vlen; ++k)
					r_mat(i, j) += xyz2(i, k) * xyz1(j, k);

		Eigen::Matrix3d sym_mat = r_mat.transpose() * r_mat;

		Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(sym_mat);
		Eigen::Vector3d eigenvalues = eigensolver.eigenvalues();
		Eigen::Matrix3d eigenvectors = eigensolver.eigenvectors();

		Eigen::Matrix3d b_mat = Eigen::Matrix3d::Zero();
		for (size_t i = 0; i < 3; ++i)
			b_mat.row(i) = r_mat * eigenvectors.col(i) / sqrt(eigenvalues[i]);

		Eigen::Matrix3d a_mat = eigenvectors.transpose();

		Eigen::Matrix3d u_mat = Eigen::Matrix3d::Zero();
		for (size_t i = 0; i < 3; ++i)
			for (size_t j = 0; j < 3; ++j)
				for (size_t k = 0; k < 3; ++k)
					u_mat(i, j) += b_mat(k, j) * a_mat(k, i);

		Eigen::Matrix3Xd fit_xyz2 = u_mat * xyz2;
		fit_xyz2.colwise() += c1;
		Eigen::Map<Eigen::VectorXd> fit_coord2(fit_xyz2.data(), fit_xyz2.size());
		return fit_coord2;
	}
	else
		return Eigen::VectorXd();
}

Eigen::VectorXd calc_displacement(Eigen::VectorXd coord1, Eigen::VectorXd coord2)
{
	if (coord1.size() == coord2.size())
		return coord2 - coord1;
	else
		return Eigen::VectorXd();
}

double calc_rmsd(Eigen::VectorXd coord1, Eigen::VectorXd coord2)
{
	if (coord1.size() == coord2.size())
	{
		Eigen::VectorXd diff = coord1 - coord2;
		size_t vlen = diff.size() / 3;

		Eigen::Map<Eigen::Array3Xd> diffxyz(diff.data(), 3, vlen);

		Eigen::VectorXd sqdist = diffxyz.pow(2).colwise().sum();
		return sqrt(sqdist.sum() / vlen);
	}
	else
		return 0.0;
}

double calc_rmsd(Eigen::VectorXd diff)
{
	size_t vlen = diff.size() / 3;

	Eigen::Map<Eigen::Array3Xd> diffxyz(diff.data(), 3, vlen);

	Eigen::VectorXd sqdist = diffxyz.pow(2).colwise().sum();
	return sqrt(sqdist.sum() / vlen);
}

double calc_mindist(Eigen::VectorXd coord1, Eigen::VectorXd coord2)
{
	Eigen::Map<Eigen::Matrix3Xd> xyz1(coord1.data(), 3, coord1.size() / 3);
	Eigen::Map<Eigen::Matrix3Xd> xyz2(coord2.data(), 3, coord2.size() / 3);
	
	Eigen::VectorXd dist(xyz1.cols());
	Eigen::Vector3d onexyz = Eigen::Vector3d::Zero();
	Eigen::Array3Xd diffxyz(3, xyz2.cols());
	for (size_t i = 0; i < size_t(xyz1.cols()); ++i)
	{
		onexyz = xyz1.col(i);
		diffxyz = xyz2.colwise() - onexyz;		
		dist(i) = diffxyz.pow(2).colwise().sum().minCoeff();
	}
	return dist.minCoeff();
}

double calc_norm(Eigen::VectorXd vector)
{
	Eigen::ArrayXd v = vector;
	return sqrt(v.pow(2).sum());
}

double calc_average_force(Eigen::VectorXd force)
{
	Eigen::Map<Eigen::Array3Xd> forcexyz(force.data(), 3, force.size() / 3);
	Eigen::VectorXd sqdist = forcexyz.pow(2).colwise().sum().sqrt();
	return sqdist.sum() / sqdist.size();
}

Eigen::MatrixXd gen_distmat(Eigen::VectorXd coord)
{
	size_t resn = coord.size() / 3;
	Eigen::MatrixXd distmat = Eigen::MatrixXd::Zero(resn, resn);
	for (size_t i = 0; i < resn; ++i)
		for (size_t j = i + 1; j < resn; ++j)
			distmat(j, i) = distmat(i, j) = sqrt(pow(coord(3 * i) - coord(3 * j), 2) + pow(coord(3 * i + 1) - coord(3 * j + 1), 2) + pow(coord(3 * i + 2) - coord(3 * j + 2), 2));

	return distmat;
}

std::list<size_t> gen_pocket(double cutoff, Eigen::VectorXd dist2ligand)
{
	std::list<size_t> pocket;
	if (cutoff > dist2ligand.minCoeff())
	{
		for (size_t i = 0; i < size_t(dist2ligand.size()); ++i)
			if (dist2ligand(i) < cutoff)
				pocket.push_back(i);
	}
	else
	{
		std::cout << std::fixed << std::setprecision(2);
		std::cout << "[Error] Given cutoff is too short. Minimum possible cutoff is " << dist2ligand.minCoeff() << "." << std::endl;
	}
	return pocket;
}

void normal_equation(Eigen::VectorXd &coeff, Eigen::MatrixXd X, Eigen::VectorXd Y)
{
	coeff = (X.transpose() * X).inverse() * X.transpose() * Y;
}

// Batch Gradient Descent
void BGD(Eigen::VectorXd &coeff, Eigen::MatrixXd X, Eigen::MatrixXd Y, double learning_rate, double convergence, size_t iterations)
{
	size_t nfeature = coeff.size();
	size_t nsample = Y.size();
	size_t nconverge = 0;
	Eigen::VectorXd tmpcoeff = Eigen::VectorXd::Zero(nfeature);
	for (size_t k = 0; k < iterations; ++k)
	{
		tmpcoeff = learning_rate / nsample * ((X * coeff - Y).transpose() * X);
		coeff -= tmpcoeff;
		for (size_t j = 0; j < nfeature; ++j)
			if (abs(tmpcoeff(j)) < convergence)
				++nconverge;
		if (nconverge == nfeature)
			break;
	}
}
