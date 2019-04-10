#pragma once
#include <list>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/QR>
#include <boost/format.hpp>

#include "manageIO.h"

using namespace std;
using namespace Eigen;

VectorXd rotate(VectorXd coord, Vector3d axis, double angle);

VectorXd fitting(VectorXd coord1, VectorXd coord2);

VectorXd calc_displacement(VectorXd coord1, VectorXd coord2);

double calc_rmsd(VectorXd coord1, VectorXd coord2);

double calc_rmsd(VectorXd diff);

double calc_mindist(VectorXd coord1, VectorXd coord2);

double calc_norm(VectorXd vector);

double calc_average_force(VectorXd force);

MatrixXd gen_distmat(VectorXd coord);

list<size_t> gen_pocket(double cutoff, VectorXd dist2ligand);

void normal_equation(VectorXd & coefficient, MatrixXd X, VectorXd Y);

void BGD(VectorXd & coeff, MatrixXd X, MatrixXd Y, double learning_rate, double convergence, size_t iterations);
