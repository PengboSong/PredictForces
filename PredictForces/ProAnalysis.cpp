#include "ProAnalysis.h"
#include "method.h"

ProAnalysis::ProAnalysis()
{
}

ProAnalysis::ProAnalysis(Pro apo, Pro binding)
{
	ProE = apo;
	ProS = binding;

	if (!ProE.empty())
	{
		hessian = ProE.gen_hessian(); // unit: J / mol
		std::cout << "[Info] Finish constructing Hessian matrix." << std::endl;
		covariance = ProE.gen_covariance(hessian); // unit: A^2 / mol
		std::cout << "[Info] Finish constructing Covariance matrix." << std::endl;

		if (!ProS.empty())
		{
			S_dist2ligand = ProS.get_dist2ligand();

			S_fitprocoord = fitting(ProE.get_procoord(), ProS.get_procoord());
			ES_displacement = calc_displacement(ProE.get_procoord(), S_fitprocoord);
			ES_force = hessian * ES_displacement;
			ES_average_force = calc_average_force(ES_force);
			ES_rmsd = calc_rmsd(ES_displacement);

			ES_info = true;
		}
	}
}

ProAnalysis::ProAnalysis(Pro apo, Pro binding, Pro allostery, Pro complex)
{
	ProE = apo;
	ProS = binding;
	ProA = allostery;
	ProAS = complex;

	if (!ProE.empty())
	{
		hessian = ProE.gen_hessian(); // unit: J / mol
		std::cout << "[Info] Finish constructing Hessian matrix." << std::endl;
		covariance = ProE.gen_covariance(hessian); // unit: A^2 / mol
		std::cout << "[Info] Finish constructing Covariance matrix." << std::endl;

		if (!ProS.empty())
		{
			S_dist2ligand = ProS.get_dist2ligand();

			S_fitprocoord = fitting(ProE.get_procoord(), ProS.get_procoord());
			ES_displacement = calc_displacement(ProE.get_procoord(), S_fitprocoord);
			ES_force = hessian * ES_displacement; // unit: J A / mol
			ES_average_force = calc_average_force(ES_force);
			ES_rmsd = calc_rmsd(ES_displacement); // unit: A

			ES_info = true;
		}
		if (!ProA.empty())
		{
			A_dist2ligand = ProA.get_dist2ligand();

			A_fitprocoord = fitting(ProE.get_procoord(), ProA.get_procoord());
			EA_displacement = calc_displacement(ProE.get_procoord(), A_fitprocoord);
			EA_force = hessian * EA_displacement; // unit: J A / mol
			EA_average_force = calc_average_force(EA_force);
			EA_rmsd = calc_rmsd(EA_displacement); // unit: A

			EA_info = true;
		}
		if (!ProAS.empty())
		{
			AS_dist2ligand = ProAS.get_dist2ligand();

			AS_fitprocoord = fitting(ProE.get_procoord(), ProAS.get_procoord());
			EAS_displacement = calc_displacement(ProE.get_procoord(), AS_fitprocoord);
			EAS_force = hessian * EAS_displacement; // unit: J A / mol
			EAS_average_force = calc_average_force(EAS_force);
			EAS_rmsd = calc_rmsd(EAS_displacement); // unit: A

			EAS_info = true;
		}
	}
}

ProAnalysis::~ProAnalysis()
{
}

Eigen::MatrixXd ProAnalysis::get_hessian()
{
	return hessian;
}

Eigen::Matrix3d ProAnalysis::get_hessian(size_t i, size_t j)
{
	if (i < size_t(hessian.rows() / 3) && j < size_t(hessian.cols() / 3))
		return Eigen::Matrix3d(hessian.block(3 * i, 3 * j, 3, 3));
	else
		return Eigen::Matrix3d();
}

double ProAnalysis::get_hessian_s(size_t si, size_t sj)
{
	if (si < size_t(hessian.rows()) && sj < size_t(hessian.cols()))
		return hessian(si, sj);
	else
		return 0.0;
}

void ProAnalysis::write_hessian(std::string writepath)
{
	std::ofstream hessianf(writepath);
	if (hessianf.is_open())
	{
		hessianf << hessian.format(CleanFmt);
		hessianf.close();
		std::cout << "[Info] Hessian matrix has been written to " << writepath << ". " << std::endl;
	}
}

Eigen::MatrixXd ProAnalysis::get_covariance()
{
	return covariance;
}

Eigen::Matrix3d ProAnalysis::get_covariance(size_t i, size_t j)
{
	if (i < size_t(covariance.rows() / 3) && j < size_t(covariance.cols() / 3))
		return Eigen::Matrix3d(covariance.block(3 * i, 3 * j, 3, 3));
	else
		return Eigen::Matrix3d();
}

double ProAnalysis::get_covariance_s(size_t si, size_t sj)
{
	if (si < size_t(covariance.rows()) && sj < size_t(covariance.cols()))
		return covariance(si, sj);
	else
		return 0.0;
}

void ProAnalysis::write_covariance(std::string writepath)
{
	std::ofstream covariancef(writepath);
	if (covariancef.is_open())
	{
		covariancef << covariance.format(CleanFmt);
		covariancef.close();
		std::cout << "[Info] Covariance matrix has been written to " << writepath << ". " << std::endl;
	}
}

void ProAnalysis::gen_free_energy()
{
	if (ES_info && EA_info)
	{
		gen_pocket_force(pocketS_force, pocketS, ES_force, ES_displacement);
		calc_energy_known(S_proenergy, S_pocketenergy, S_energy, pocketS_force, ProS.get_distmat());

		gen_pocket_force(pocketA_force, pocketA, EA_force, EA_displacement);
		calc_energy_known(A_proenergy, A_pocketenergy, A_energy, pocketA_force, ProA.get_distmat());

		pocketAS_force = pocketS_force + pocketA_force; // Right?
		calc_energy_unknown(AS_proenergy, AS_pocketenergy, AS_energy, pocketAS_force);
		ddG = AS_energy - S_energy - A_energy;

		print_energy_results();
	}
	else if (EAS_info && ES_info)
	{
		gen_pocket_force(pocketS_force, pocketS, ES_force, ES_displacement);
		calc_energy_known(S_proenergy, S_pocketenergy, S_energy, pocketS_force, ProS.get_distmat());

		gen_pocket_force(pocketAS_force, pocketAS, EAS_force, EAS_displacement);
		calc_energy_known(AS_proenergy, AS_pocketenergy, AS_energy, pocketAS_force, ProAS.get_distmat());

		pocketA_force = Eigen::VectorXd::Zero(pocketAS_force.size());
		for (std::list<size_t>::iterator it = pocketAS.begin(); it != pocketAS.end(); ++it)
		{
			if (!in_pocketS(*it))
			{
				add_to_pocketA(*it);
				pocketA_force(*it * 3) = pocketAS_force(*it * 3);
				pocketA_force(*it * 3 + 1) = pocketAS_force(*it * 3 + 1);
				pocketA_force(*it * 3 + 2) = pocketAS_force(*it * 3 + 2);
			}
		}

		calc_energy_unknown(A_proenergy, A_pocketenergy, A_energy, pocketA_force);
		ddG = AS_energy - S_energy - A_energy;

		print_energy_results();
	}
	else if (EAS_info && EA_info)
	{
		gen_pocket_force(pocketA_force, pocketA, EA_force, EA_displacement);
		calc_energy_known(A_proenergy, A_pocketenergy, A_energy, pocketA_force, ProA.get_distmat());

		gen_pocket_force(pocketAS_force, pocketAS, EAS_force, EAS_displacement);
		calc_energy_known(AS_proenergy, AS_pocketenergy, AS_energy, pocketAS_force, ProAS.get_distmat());

		pocketS_force = Eigen::VectorXd::Zero(pocketAS_force.size());
		for (std::list<size_t>::iterator it = pocketAS.begin(); it != pocketAS.end(); ++it)
		{
			if (!in_pocketA(*it))
			{
				add_to_pocketS(*it);
				pocketS_force(*it * 3) = pocketAS_force(*it * 3);
				pocketS_force(*it * 3 + 1) = pocketAS_force(*it * 3 + 1);
				pocketS_force(*it * 3 + 2) = pocketAS_force(*it * 3 + 2);
			}
		}

		calc_energy_unknown(S_proenergy, S_pocketenergy, S_energy, pocketS_force);
		ddG = AS_energy - S_energy - A_energy;

		print_energy_results();
	}
	else
		std::cout << "[Error] Lack necessary information." << std::endl;
}

double ProAnalysis::calc_model_rmsd(std::list<size_t> pocket, Eigen::VectorXd pocket_force, Eigen::VectorXd pro_force, Eigen::VectorXd refcoord, Eigen::VectorXd displacenment)
{
	gen_pocket_force(pocket_force, pocket, pro_force, displacenment);

	Eigen::VectorXd mprocoord = covariance * pocket_force + ProE.get_procoord();
	Eigen::VectorXd fitmprocoord = fitting(refcoord, mprocoord);
	return calc_rmsd(refcoord, fitmprocoord);
}

void ProAnalysis::show_pocket(std::list<size_t> pocket)
{
	std::cout << "[Info] Pocket residues: ";
	for (std::list<size_t>::iterator it = pocket.begin(); it != pocket.end(); ++it)
	{
		std::cout << *it << " ";
	}
	std::cout << std::endl;
}

void ProAnalysis::test_pocket(bool info, std::list<size_t> pocket, Eigen::VectorXd pocket_force, Eigen::VectorXd pro_force, Eigen::VectorXd refcoord, Eigen::VectorXd displacenment)
{
	if (!pocket.empty())
	{
		show_pocket(pocket);

		if (info)
		{
			std::cout << std::fixed << std::setprecision(4);
			std::cout << "[Info] RMSD between real structure and structure calculated according to current pocket: " << calc_model_rmsd(pocket, pocket_force, pro_force, refcoord, displacenment) << " A." << std::endl;
		}
	}
	else
		std::cout << "[Error] The binding pocket domain is not specificed." << std::endl;
}

void ProAnalysis::show_pocket_force(std::list<size_t> pocket, Eigen::VectorXd pocket_force)
{
	for (std::list<size_t>::iterator it = pocket.begin(); it != pocket.end(); ++it)
	{
		Eigen::Vector3d resforce = Eigen::Vector3d::Zero();
		resforce << pocket_force(*it * 3), pocket_force(*it * 3 + 1), pocket_force(*it * 3 + 2);
		std::cout << "RES " << *it << " FORCE " << calc_norm(resforce) << std::endl;
	}
}

bool ProAnalysis::in_pocket(std::list<size_t> pocket, size_t id)
{
	bool find_element_flag = false;
	for (std::list<size_t>::iterator it = pocket.begin(); it != pocket.end(); ++it)
		if (*it == id)
			find_element_flag = true;
	return find_element_flag;
}

void ProAnalysis::add_to_pocket(std::list<size_t> pocket, size_t id)
{
	if (id < ProE.get_resn())
	{
		if (in_pocket(pocket, id))
			std::cout << "[Info] Given residue ID already in pocket." << std::endl;
		else
			pocket.push_back(id);
	}
	else
		std::cout << "[Error] Given residue ID out of range." << std::endl;
}

void ProAnalysis::remove_from_pocket(std::list<size_t> pocket, size_t id)
{
	if (id < ProE.get_resn())
	{
		if (in_pocket(pocket, id))
		{
			for (std::list<size_t>::iterator it = pocket.begin(); it != pocket.end(); ++it)
				if (*it == id)
					pocket.erase(it);
		}
		else
			std::cout << "[Info] Given residue ID not in pocket. Can not erase it." << std::endl;
	}
	else
		std::cout << "[Error] Given residue ID out of range." << std::endl;
}

void ProAnalysis::gen_pocket(bool has_ligand, std::list<size_t> &pocket, double cutoff, Eigen::VectorXd dist2ligand)
{
	if (has_ligand)
	{
		if (cutoff > dist2ligand.minCoeff())
		{
			if (!pocket.empty())
				pocket.clear();
			for (size_t i = 0; i < size_t(dist2ligand.size()); ++i)
				if (dist2ligand(i) < cutoff)
					pocket.push_back(i);
		}
		else
		{
			std::cout << std::fixed << std::setprecision(2);
			std::cout << "[Error] Given cutoff is too short. Minimum possible cutoff is " << dist2ligand.minCoeff() << "." << std::endl;
		}
	}
	else
		std::cout << "[Error] Can not find ligand information." << std::endl;
}

void ProAnalysis::gen_pocket_force(Eigen::VectorXd &pocket_force, std::list<size_t> pocket, Eigen::VectorXd pro_force, Eigen::VectorXd displacenment)
{
	pocket_force = Eigen::VectorXd::Zero(covariance.rows());

	size_t ndim = pocket.size() * 3;
	Eigen::MatrixXd X = Eigen::MatrixXd::Zero(covariance.rows(), ndim);
	Eigen::VectorXd Y = displacenment;
	Eigen::VectorXd coeff = Eigen::VectorXd::Zero(ndim);

	size_t i = 0;
	for (std::list<size_t>::iterator it = pocket.begin(); it != pocket.end(); ++it)
	{
		X.col(3 * i) = covariance.col(*it * 3);
		X.col(3 * i + 1) = covariance.col(*it * 3 + 1);
		X.col(3 * i + 2) = covariance.col(*it * 3 + 2);

		coeff(i * 3) = pro_force(*it * 3);
		coeff(i * 3 + 1) = pro_force(*it * 3 + 1);
		coeff(i * 3 + 2) = pro_force(*it * 3 + 2);

		++i;
	}

	/*
	std::cout << std::resetiosflags(std::ios::fixed);
	std::cout << "Using learning step " << LEARNING_STEP << "." << std::endl;
	std::cout << "Using convergence " << CONVERGENCE << "." << std::endl;
	std::cout << std::setiosflags(std::ios::fixed);
	BGD(coeff, X, Y, LEARNING_STEP, CONVERGENCE, ITERATION_TIMES);
	*/
	normal_equation(coeff, X, Y);

	i = 0;
	for (std::list<size_t>::iterator it = pocket.begin(); it != pocket.end(); ++it)
	{
		pocket_force(*it * 3) = coeff(i * 3);
		pocket_force(*it * 3 + 1) = coeff(i * 3 + 1);
		pocket_force(*it * 3 + 2) = coeff(i * 3 + 2);
	}
	
}

void ProAnalysis::calc_energy_known(double &proenergy, double &pocketenergy, double &totenergy, Eigen::VectorXd pocket_force, Eigen::MatrixXd distmat)
{
	Eigen::ArrayXXd distdiffmat = distmat - ProE.get_distmat();
	proenergy = (distdiffmat.pow(2) * ProE.get_kmat()).sum();
	pocketenergy = -pocket_force.transpose() * covariance * pocket_force;
	totenergy = proenergy + pocketenergy;
}

void ProAnalysis::calc_energy_unknown(double &proenergy, double &pocketenergy, double &totenergy, Eigen::VectorXd pocket_force)
{
	Eigen::VectorXd procoord = covariance * pocket_force + ProE.get_procoord();
	Eigen::ArrayXXd distdiffmat = gen_distmat(procoord) - ProE.get_distmat();
	proenergy = (distdiffmat.pow(2) * ProE.get_kmat()).sum();
	pocketenergy = -pocket_force.transpose() * covariance * pocket_force;
	totenergy = proenergy + pocketenergy;
}

void ProAnalysis::debug_energy_unknown(Eigen::VectorXd pocket_force)
{
	Eigen::VectorXd procoord = covariance * pocket_force + ProE.get_procoord();
	Eigen::ArrayXXd distdiffmat = gen_distmat(procoord) - ProE.get_distmat();
	double proenergy = (distdiffmat.pow(2) * ProE.get_kmat()).sum();
	double pocketenergy = -pocket_force.transpose() * covariance * pocket_force;
	double totenergy = proenergy + pocketenergy;

	std::cout << "---< DEBUG >---" << std::endl;
	std::cout << std::fixed << std::setprecision(4);
	std::cout << "Total: " << totenergy * 1e-3 << std::endl;
	std::cout << "Pro: " << proenergy * 1e-3 << " kJ/mol." << std::endl;
	std::cout << "Pocket: " << proenergy * 1e-3 << " kJ/mol." << std::endl;
	std::cout << "---< END >---" << std::endl;
}

void ProAnalysis::print_energy_results()
{
	std::cout << "[Info] Free energy results: " << std::endl;
	std::cout << std::fixed << std::setprecision(4);
	std::cout << "Free energy for binding state structure S: " << S_energy * 1e-3 << " kJ/mol." << std::endl;
	std::cout << "Pro: " << S_proenergy * 1e-3 << " kJ/mol." << std::endl;
	std::cout << "Pocket: " << S_pocketenergy * 1e-3 << " kJ/mol." << std::endl;
	std::cout << "Free energy for allostery state structure A: " << A_energy * 1e-3 << " kJ/mol." << std::endl;
	std::cout << "Pro: " << A_proenergy * 1e-3 << " kJ/mol." << std::endl;
	std::cout << "Pocket: " << A_pocketenergy * 1e-3 << " kJ/mol." << std::endl;
	std::cout << "Free energy for complex structure AS: " << AS_energy * 1e-3 << " kJ/mol." << std::endl;
	std::cout << "Pro: " << AS_proenergy * 1e-3 << " kJ/mol." << std::endl;
	std::cout << "Pocket: " << AS_pocketenergy * 1e-3 << " kJ/mol." << std::endl;
	std::cout << "Change of free energy: " << ddG * 1e-3 << " kJ/mol." << std::endl;

}

