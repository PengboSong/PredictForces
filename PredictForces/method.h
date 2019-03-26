#pragma once

#include <iostream>

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
