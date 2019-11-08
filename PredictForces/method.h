#pragma once
#include <list>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/QR>
#include <boost/format.hpp>

#include "defines.h"
#include "handle_io.h"

using namespace Eigen;

Matrix3Xd coord_vec2mat(VectorXd coord);

VectorXd coord_mat2vec(Matrix3Xd coord_xyz);

VectorXd rotate(VectorXd coord, Vector3d axis, double angle);

VectorXd fitting(VectorXd coord1, VectorXd coord2);

VectorXd calc_displacement(VectorXd coord1, VectorXd coord2);

Vector3d coord_center(VectorXd coord);

double calc_rmsd(VectorXd coord1, VectorXd coord2);

double calc_rmsd(VectorXd diff);

double calc_mindist(VectorXd coord1, VectorXd coord2);

double calc_norm(VectorXd vector);

double calc_average(VectorXd force);

MatrixXd gen_distmat(DistMatType t, VectorXd coord, double cutoff = 0.0);

MatrixXd gen_differ(MatrixXd X, MatrixXd Y);

std::list<size_t> gen_pocket(double cutoff, VectorXd dist2ligand);

Matrix3d euler_rotation_matrix(Vector3d degrees);

void grep_pocket_coord(VectorXd &pocket_coord, VectorXd original_coord, PocketList pocket_members);

void modify_pocket_coord(VectorXd &coord, VectorXd replace_coord, PocketList pocket_members);

void copy_pocket_coord(VectorXd &coord, VectorXd source_coord, PocketList pocket_members);

void normal_equation(VectorXd & coefficient, MatrixXd X, VectorXd Y);

void BGD(VectorXd & coeff, MatrixXd X, VectorXd Y, BGDpara paras, bool checkinf = false);
