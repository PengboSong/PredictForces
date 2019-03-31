#pragma once
#include <iostream>
#include <iomanip>
#include <list>
#include <cmath>

#include <Eigen/Dense>
#include <Eigen/QR>

Eigen::VectorXd rotate(Eigen::VectorXd coord, Eigen::Vector3d axis, double angle);

Eigen::VectorXd fitting(Eigen::VectorXd coord1, Eigen::VectorXd coord2);

Eigen::VectorXd calc_displacement(Eigen::VectorXd coord1, Eigen::VectorXd coord2);

double calc_rmsd(Eigen::VectorXd coord1, Eigen::VectorXd coord2);

double calc_rmsd(Eigen::VectorXd diff);

double calc_mindist(Eigen::VectorXd coord1, Eigen::VectorXd coord2);

double calc_norm(Eigen::VectorXd vector);

double calc_average_force(Eigen::VectorXd force);

Eigen::MatrixXd gen_distmat(Eigen::VectorXd coord);

std::list<size_t> gen_pocket(double cutoff, Eigen::VectorXd dist2ligand);

void normal_equation(Eigen::VectorXd & coefficient, Eigen::MatrixXd X, Eigen::VectorXd Y);

void BGD(Eigen::VectorXd & coeff, Eigen::MatrixXd X, Eigen::MatrixXd Y, double learning_rate, double convergence, size_t iterations);
