#pragma once
#include <list>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/QR>
#include <boost/format.hpp>

#include "handle_io.h"

using namespace Eigen;

enum DistMatType : uint8_t {
	Dist,
	XYZdiff,
	XYZdiffLtCutoff
};

typedef struct {
	double learning_rate = 0.0;
	double convergence = 0.0;
	size_t niteration = 0;
} BGDpara;

VectorXd rotate(VectorXd coord, Vector3d axis, double angle);

VectorXd fitting(VectorXd coord1, VectorXd coord2);

VectorXd calc_displacement(VectorXd coord1, VectorXd coord2);

double calc_rmsd(VectorXd coord1, VectorXd coord2);

double calc_rmsd(VectorXd diff);

double calc_mindist(VectorXd coord1, VectorXd coord2);

double calc_norm(VectorXd vector);

double calc_average(VectorXd force);

MatrixXd gen_distmat(DistMatType t, VectorXd coord, double cutoff = 0.0);

MatrixXd gen_differ(MatrixXd X, MatrixXd Y);

std::list<size_t> gen_pocket(double cutoff, VectorXd dist2ligand);

Matrix3d euler_rotation_matrix(Vector3d degrees);

VectorXd coords_derivate(VectorXd coord, MatrixXd distmat_0, MatrixXd distmat, ArrayXXd kmat);

Vector3d euler_degrees_derivate(Vector3d r, Vector3d dvdr, Vector3d degrees);

VectorXd degrees_derivate(VectorXd coord, VectorXd dR, Vector3d vcenter, Vector3d degrees);

void normal_equation(VectorXd & coefficient, MatrixXd X, VectorXd Y);

void BGD(VectorXd & coeff, MatrixXd X, VectorXd Y, BGDpara paras, bool checkinf = false);
