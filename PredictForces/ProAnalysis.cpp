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
		handle_info("Finish constructing Hessian matrix.");
		covariance = ProE.gen_covariance(hessian); // unit: A^2 / mol
		handle_info("Finish constructing Covariance matrix.");

		if (!ProS.empty())
		{
			S_dist2ligand = ProS.get_dist2ligand();

			S_fitprocoord = fitting(ProE.get_procoord(), ProS.get_procoord());
			ES_displacement = calc_displacement(ProE.get_procoord(), S_fitprocoord);
			ES_force = hessian * ES_displacement;
			// ES_average_force = calc_average_force(ES_force);
			ES_rmsd = calc_rmsd(ES_displacement);
			std::cout << "rmsd from pdb: " << "\t" << ES_rmsd << std::endl;
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
		handle_info("Finish constructing Hessian matrix.");
		//write_hessian("C:\\Users\\hqj\\Desktop\\hessian");
		covariance = ProE.gen_covariance(hessian); // unit: A^2 / mol
		handle_info("Finish constructing Covariance matrix.");
		//write_covariance("C:\\Users\\hqj\\Desktop\\covariance");

		if (!ProS.empty())
		{
			if (ProE.get_resn() != ProS.get_resn())
				handle_error("Apo state protein and binding state protein do not have equal residue numbers.");

			S_dist2ligand = ProS.get_dist2ligand();

			S_fitprocoord = fitting(ProE.get_procoord(), ProS.get_procoord());
			handle_info("Fitting process succeed.");
			ES_displacement = calc_displacement(ProE.get_procoord(), S_fitprocoord);
			handle_info("Calculating displacement succeed.");
			ES_force = hessian * ES_displacement; // unit: J A / mol
			handle_info("Calculating force succeed.");
			//std::cout << "originES_force: " << "\n" << ES_force << std::endl;
			// ES_average_force = calc_average_force(ES_force);
			ES_rmsd = calc_rmsd(ES_displacement); // unit: A
			handle_result(boost::format("RMSD from PDB file: %1% A.") % EA_rmsd);
			
			ES_info = true;
		}
		if (!ProA.empty())
		{
			if (ProE.get_resn() != ProA.get_resn())
				handle_error("Apo state protein and allostery state protein do not have equal residue numbers.");

			A_dist2ligand = ProA.get_dist2ligand();

			A_fitprocoord = fitting(ProE.get_procoord(), ProA.get_procoord());
			handle_info("Fitting process succeed.");
			EA_displacement = calc_displacement(ProE.get_procoord(), A_fitprocoord);
			handle_info("Calculating displacement succeed.");
			EA_force = hessian * EA_displacement; // unit: J A / mol
			handle_info("Calculating force succeed.");
			//std::cout << "originEA_force: " << "\n" << EA_force << std::endl;
			// EA_average_force = calc_average_force(EA_force);
			EA_rmsd = calc_rmsd(EA_displacement); // unit: A
			handle_result(boost::format("RMSD from PDB file: %1% A.") % EA_rmsd);

			EA_info = true;
		}
		if (!ProAS.empty())
		{
			if (ProE.get_resn() != ProA.get_resn())
				handle_error("Apo state protein and complex state protein do not have equal residue numbers.");
			
			AS_dist2ligand = ProAS.get_dist2ligand();

			AS_fitprocoord = fitting(ProE.get_procoord(), ProAS.get_procoord());
			handle_info("Fitting process succeed.");
			EAS_displacement = calc_displacement(ProE.get_procoord(), AS_fitprocoord);
			handle_info("Calculating displacement succeed.");
			EAS_force = hessian * EAS_displacement; // unit: J A / mol
			handle_info("Calculating force succeed.");
			//std::cout << "originEAS_force: " << "\n" << EAS_force << std::endl;
			// EAS_average_force = calc_average_force(EAS_force);
			EAS_rmsd = calc_rmsd(EAS_displacement); // unit: A
			handle_result(boost::format("RMSD from PDB file: %1% A.") % EA_rmsd);

			EAS_info = true;
		}
	}
}

ProAnalysis::~ProAnalysis()
{
}

void ProAnalysis::interactive_pocket(unsigned int mode)
{
	string buf, cmd;
	vector<string> para;
	while (true)
	{
		cout << '[' << mode << ']' << " >>> ";
		getline(cin, buf);
		trim(buf);
		split(para, buf, boost::is_any_of(" "));
		cmd = para[0];
		if (cmd == "back")
			break;
		else if (cmd == "add")
		{
			if (!para[1].empty())
			{
				size_t resid = lexical_cast<size_t>(para[1]) - 1;
				switch (mode)
				{
				case 0:
					add_to_pocketS(resid);
					break;
				case 1:
					add_to_pocketA(resid);
					break;
				case 2:
					add_to_pocketAS(resid);
					break;
				}
			}
			else
				handle_hint("No residue ID given. Please enter again.");
		}
		else if (cmd == "del")
		{
			if (!para[1].empty())
			{
				size_t resid = lexical_cast<size_t>(para[1]) - 1;
				switch (mode)
				{
				case 0:
					remove_from_pocketS(resid);
					break;
				case 1:
					remove_from_pocketA(resid);
					break;
				case 2:
					remove_from_pocketAS(resid);
					break;
				}
			}
			else
				handle_hint("No residue ID given. Please enter again.");
		}
		else if (cmd == "gen-pocket")
		{
			if (!para[1].empty())
			{
				double cutoff = lexical_cast<double>(para[1]);
				switch (mode)
				{
				case 0:
					gen_pocketS(cutoff);
					break;
				case 1:
					gen_pocketA(cutoff);
					break;
				case 2:
					gen_pocketAS(cutoff);
					break;
				}
			}
			else
				handle_hint("No cutoff length given. Please enter again.");
		}
		else if (cmd == "show")
		{
			switch (mode)
			{
			case 0:
				show_pocketS();
				break;
			case 1:
				show_pocketA();
				break;
			case 2:
				show_pocketAS();
				break;
			}
		}
		else if (cmd == "test")
		{
			switch (mode)
			{
			case 0:
				test_pocketS();
				break;
			case 1:
				test_pocketA();
				break;
			case 2:
				test_pocketAS();
				break;
			}
		}
		else if (cmd == "gen-force")
		{
			switch (mode)
			{
			case 0:
				gen_pocketS_force();
				break;
			case 1:
				gen_pocketA_force();
				break;
			case 2:
				gen_pocketAS_force();
				break;
			}
		}
		else if (cmd == "show-force")
		{
			switch (mode)
			{
			case 0:
				show_pocketS_force();
				break;
			case 1:
				show_pocketA_force();
				break;
			case 2:
				show_pocketAS_force();
				break;
			}
		}
		else if (cmd == "origin-force")
		{
			switch (mode)
			{
			case 0:
				show_pro_pocketS_force();
				break;
			case 1:
				show_pro_pocketA_force();
				break;
			case 2:
				show_pro_pocketAS_force();
				break;
			}
		}
		else if (cmd == "all-origin-force")
		{
			switch (mode)
			{
			case 0:
				show_proS_all_force();
				break;
			case 1:
				show_proA_all_force();
				break;
			case 2:
				show_proAS_all_force();
				break;
			}
		}
		else
			handle_hint("Unknown command.");
		cmd.clear();
	}
}

void ProAnalysis::interactive()
{
	string buf, cmd;
	vector<string> para;
	while (true)
	{
		cout << ">>> ";
		getline(cin, buf);
		trim(buf);
		split(para, buf, boost::is_any_of(" "));
		cmd = para[0];
		if (cmd == "exit")
			break;
		else if (cmd == "pocketS")
			interactive_pocket(0);
		else if (cmd == "pocketA")
			interactive_pocket(1);
		else if (cmd == "pocketAS")
			interactive_pocket(2);
		else if (cmd == "energy")
			gen_free_energy();
		else if (cmd == "LFmethod")
		{
			choose_LFmethod();
			show_LFmethod_detail();
		}
		else if (cmd == "LFpara")
		{
			// TODO;
		}
		else
			handle_hint("Unknown command.");
		cmd.clear();
	}
}

MatrixXd ProAnalysis::get_hessian()
{
	return hessian;
}

Matrix3d ProAnalysis::get_hessian(size_t i, size_t j)
{
	if (i < size_t(hessian.rows() / 3) && j < size_t(hessian.cols() / 3))
		return Matrix3d(hessian.block(3 * i, 3 * j, 3, 3));
	else
		return Matrix3d();
}

double ProAnalysis::get_hessian_s(size_t si, size_t sj)
{
	if (si < size_t(hessian.rows()) && sj < size_t(hessian.cols()))
		return hessian(si, sj);
	else
		return 0.0;
}


void ProAnalysis::write_hessian(string writepath)
{
	ofstream hessianf(writepath);
	if (hessianf.is_open())
	{
		hessianf << hessian.format(CleanFmt);
		hessianf.close();
		handle_info(boost::format("Hessian matrix has been written to %1%.") % writepath);
	}
}


MatrixXd ProAnalysis::get_covariance()
{
	return covariance;
}

Matrix3d ProAnalysis::get_covariance(size_t i, size_t j)
{
	if (i < size_t(covariance.rows() / 3) && j < size_t(covariance.cols() / 3))
		return Matrix3d(covariance.block(3 * i, 3 * j, 3, 3));
	else
		return Matrix3d();
}

double ProAnalysis::get_covariance_s(size_t si, size_t sj)
{
	if (si < size_t(covariance.rows()) && sj < size_t(covariance.cols()))
		return covariance(si, sj);
	else
		return 0.0;
}

void ProAnalysis::write_covariance(string writepath)
{
	ofstream covariancef(writepath);
	if (covariancef.is_open())
	{
		covariancef << covariance.format(CleanFmt);
		covariancef.close();
		handle_info(boost::format("Covariance matrix has been written to %1%.") % writepath);
	}
}

void ProAnalysis::show_LFmethod_detail()
{
	switch (LFmethod_mode)
	{
	case 0:
		// Normal Equation
		handle_info("Multiple Linear Fitting using normal equation.");
		break;
	case 1:
		// Batch Gradient Descent
		handle_info("Multiple Linear Fitting using batch gradient descent.");
		handle_info(boost::format("Using learning step %1$.2e") % LEARNING_STEP);
		handle_info(boost::format("Using convergence %1$.2e") % CONVERGENCE);
		handle_info(boost::format("Using maximum iteration times %1%") % ITERATION_TIMES);
		handle_info(boost::format("Using random times %1%") % RANDOM_TIMES);
		break;
	default:
		handle_warning("Invalid LF method mode.");
	}
}

void ProAnalysis::set_LFmethod(unsigned int mode)
{
	if (LFmethods.find(mode) != LFmethods.end())
		LFmethod_mode = mode;
}

void ProAnalysis::choose_LFmethod()
{
	for (map<unsigned int, string>::iterator it = LFmethods.begin(); it != LFmethods.end(); ++it)
		handle_hint(boost::format("%2% - %1%") % it->first % it->second);
	string buf;
	cout << "Enter method id:";
	cin >> buf;
	set_LFmethod(lexical_cast<unsigned int>(buf));
}

void ProAnalysis::gen_free_energy()
{
	if (ES_info && EA_info)
	{
		//gen_pocket_force(has_pocketS_force_flag, pocketS_force, pocketS, ES_force, ES_displacement);
		calc_energy_known(has_pocketS_force_flag, S_proenergy, S_pocketenergy, S_energy, pocketS, ES_force, ProS.get_distmat(), ES_displacement);
		calc_energy_unknown(has_pocketS_force_flag, S_predict_proenergy, S_predict_pocketenergy, S_predict_energy, pocketS_force);

		//gen_pocket_force(has_pocketA_force_flag, pocketA_force, pocketA, EA_force, EA_displacement);
		calc_energy_known(has_pocketA_force_flag, A_proenergy, A_pocketenergy, A_energy, pocketA, EA_force, ProA.get_distmat(), EA_displacement);
		calc_energy_unknown(has_pocketA_force_flag, A_predict_proenergy, A_predict_pocketenergy, A_predict_energy, pocketA_force);

		pocketAS_force = pocketS_force + pocketA_force;
		has_pocketAS_force_flag = true;
		calc_energy_unknown(has_pocketAS_force_flag, AS_proenergy, AS_pocketenergy, AS_energy, pocketAS_force);
		AS_predict_proenergy = AS_proenergy;
		AS_predict_pocketenergy = AS_pocketenergy;
		AS_predict_energy = AS_energy;
		ddG = AS_energy - S_energy - A_energy;
		ddG_predict = AS_predict_energy - S_predict_energy - A_predict_energy;

		print_energy_results();
	}
	else if (EAS_info && ES_info)
	{
		//gen_pocket_force(has_pocketS_force_flag, pocketS_force, pocketS, ES_force, ES_displacement);
		calc_energy_known(has_pocketS_force_flag, S_proenergy, S_pocketenergy, S_energy, pocketS, ES_force, ProS.get_distmat(), ES_displacement);
		calc_energy_unknown(has_pocketS_force_flag, S_predict_proenergy, S_predict_pocketenergy, S_predict_energy, pocketS_force);

		// gen_pocket_force(has_pocketAS_force_flag, pocketAS_force, pocketAS, EAS_force, EAS_displacement);
		//gen_pocket_force(has_pocketAS_force_flag, pocketAS_force, pocketS_force, pocketAS, pocketS, EAS_force, EAS_displacement);
		calc_energy_known(has_pocketAS_force_flag, AS_proenergy, AS_pocketenergy, AS_energy, pocketAS, EAS_force, ProAS.get_distmat(), EAS_displacement);
		calc_energy_unknown(has_pocketAS_force_flag, AS_predict_proenergy, AS_predict_pocketenergy, AS_predict_energy, pocketAS_force);

		pocketA_force = VectorXd::Zero(pocketAS_force.size());
		for (list<size_t>::iterator it = pocketAS.begin(); it != pocketAS.end(); ++it)
		{
			if (!in_pocketS(*it))
			{
				add_to_pocketA(*it);
				pocketA_force(*it * 3) = pocketAS_force(*it * 3);
				pocketA_force(*it * 3 + 1) = pocketAS_force(*it * 3 + 1);
				pocketA_force(*it * 3 + 2) = pocketAS_force(*it * 3 + 2);
			}
		}
		has_pocketA_force_flag = true;

		calc_energy_unknown(has_pocketA_force_flag, A_proenergy, A_pocketenergy, A_energy, pocketA_force);
		A_predict_proenergy = A_proenergy;
		A_predict_pocketenergy = A_pocketenergy;
		A_predict_energy = A_energy;
		ddG = AS_energy - S_energy - A_energy;
		ddG_predict = AS_predict_energy - S_predict_energy - A_predict_energy;

		print_energy_results();
	}
	else if (EAS_info && EA_info)
	{
		//gen_pocket_force(has_pocketA_force_flag, pocketA_force, pocketA, EA_force, EA_displacement);
		calc_energy_known(has_pocketA_force_flag, A_proenergy, A_pocketenergy, A_energy, pocketA, EA_force, ProA.get_distmat(), EA_displacement);
		calc_energy_unknown(has_pocketA_force_flag, A_predict_proenergy, A_predict_pocketenergy, A_predict_energy, pocketA_force);

		// gen_pocket_force(has_pocketAS_force_flag, pocketAS_force, pocketAS, EAS_force, EAS_displacement);
		//gen_pocket_force(has_pocketAS_force_flag, pocketAS_force, pocketA_force, pocketAS, pocketA, EAS_force, EAS_displacement);
		calc_energy_known(has_pocketAS_force_flag, AS_proenergy, AS_pocketenergy, AS_energy, pocketAS, EAS_force, ProAS.get_distmat(), EAS_displacement);
		calc_energy_unknown(has_pocketAS_force_flag, AS_predict_proenergy, AS_predict_pocketenergy, AS_predict_energy, pocketAS_force);

		pocketS_force = VectorXd::Zero(pocketAS_force.size());
		for (list<size_t>::iterator it = pocketAS.begin(); it != pocketAS.end(); ++it)
		{
			if (!in_pocketA(*it))
			{
				add_to_pocketS(*it);
				pocketS_force(*it * 3) = pocketAS_force(*it * 3);
				pocketS_force(*it * 3 + 1) = pocketAS_force(*it * 3 + 1);
				pocketS_force(*it * 3 + 2) = pocketAS_force(*it * 3 + 2);
			}
		}
		has_pocketS_force_flag = true;

		calc_energy_unknown(has_pocketS_force_flag, S_proenergy, S_pocketenergy, S_energy, pocketS_force);
		S_predict_proenergy = S_proenergy;
		S_predict_pocketenergy = S_pocketenergy;
		S_predict_energy = S_energy;
		ddG = AS_energy - S_energy - A_energy;
		ddG_predict = AS_predict_energy - S_predict_energy - A_predict_energy;

		print_energy_results();
	}
	else
		handle_error("Lack necessary protein information.");
}

void ProAnalysis::switch_LFmethod(VectorXd & coeff, MatrixXd X, VectorXd Y)
{
	switch (LFmethod_mode)
	{
	case 0:
		normal_equation(coeff, X, Y);
	case 1:
		BGD(coeff, X, Y, LEARNING_STEP, CONVERGENCE, ITERATION_TIMES, RANDOM_TIMES);
	default:
		BGD(coeff, X, Y, LEARNING_STEP, CONVERGENCE, ITERATION_TIMES, RANDOM_TIMES);
	}
}

double ProAnalysis::calc_model_rmsd(bool flag, VectorXd pocket_force, VectorXd refcoord)
{
	if (flag)
	{
		VectorXd mprocoord = covariance * pocket_force / kB / Temp / Navo + ProE.get_procoord();
		VectorXd fitmprocoord = fitting(refcoord, mprocoord);
		return calc_rmsd(refcoord, fitmprocoord);
	}
	else
	{
		handle_warning("Can not calculate model RMSD without pocket force generated. Call \"gen_pocket*_force\" function first.");
		return 0.0;
	}
}

void ProAnalysis::show_pocket(list<size_t> pocket)
{
	string buf = "Pocket residues: ";
	for (list<size_t>::iterator it = pocket.begin(); it != pocket.end(); ++it)
	{
		buf += to_string(*it + 1) + " ";
	}

	handle_result(buf);
}

void ProAnalysis::test_pocket(bool flag, bool info, list<size_t> pocket, VectorXd pocket_force, VectorXd refcoord)
{
	if (!pocket.empty() && flag)
	{
		show_pocket(pocket);
		if (info)
			handle_result(boost::format("RMSD between real structure and structure calculated according to current pocket: %1% A.") % calc_model_rmsd(flag, pocket_force, refcoord));
		else
			handle_warning("Lack necessary protein information.");
	}
	else
		handle_warning("The binding pocket domain is not specificed.");
}

void ProAnalysis::show_pocket_force(bool flag, list<size_t> pocket, VectorXd pocket_force)
{
	vector<string> buf;
	if (flag)
	{
		for (list<size_t>::iterator it = pocket.begin(); it != pocket.end(); ++it)
		{
			Vector3d resforce = Vector3d::Zero();
			resforce << pocket_force(*it * 3), pocket_force(*it * 3 + 1), pocket_force(*it * 3 + 2);
			buf.push_back(
				(boost::format("RES %1% FORCE %2$.4f") % (*it + 1) % calc_norm(resforce)).str()
			);
		}
		handle_result("Pocket force:", buf);
	}
	else
		handle_warning("The binding pocket domain is not specificed.");
}

void ProAnalysis::show_pro_pocket_force(list<size_t> pocket, VectorXd pro_force)
{
	vector<string> buf;
	for (list<size_t>::iterator it = pocket.begin(); it != pocket.end(); ++it)
	{
		Vector3d resforce = Vector3d::Zero();
		resforce << pro_force(*it * 3), pro_force(*it * 3 + 1), pro_force(*it * 3 + 2);
		buf.push_back(
			(boost::format("RES %1% FORCE %2$.4f") % (*it + 1) % calc_norm(resforce)).str()
		);
	}
	handle_result("Pocket force (from structure):", buf);
}

void ProAnalysis::show_pro_all_force(VectorXd pro_force)
{
	vector<string> buf;
	vector<pair<size_t, double>> forces;
	for (size_t i = 0; i < ProE.get_resn(); ++i)
	{
		Vector3d resforce = Vector3d::Zero();
		resforce << pro_force(i * 3), pro_force(i * 3 + 1), pro_force(i * 3 + 2);
		forces.push_back(make_pair(i, calc_norm(resforce)));
	}
	sort(forces.begin(), forces.end(), [](pair<size_t, double> a, pair<size_t, double> b) {
		return a.second > b.second;
	});

	for (size_t i = 0; i < ProE.get_resn(); ++i)
		buf.push_back(
			(boost::format("RES %1% FORCE %2$.4f") % (forces[i].first + 1) % forces[i].second).str()
		);

	handle_result("Protein force", buf);
}

bool ProAnalysis::in_pocket(list<size_t> pocket, size_t id)
{
	bool find_element_flag = false;
	for (list<size_t>::iterator it = pocket.begin(); it != pocket.end(); ++it)
		if (*it == id)
			find_element_flag = true;
	return find_element_flag;
}

void ProAnalysis::add_to_pocket(list<size_t> & pocket, size_t id)
{
	if (id < ProE.get_resn())
	{
		if (in_pocket(pocket, id))
			handle_info("Given residue ID already in pocket.");
		else
		{
			pocket.push_back(id);
			handle_info(boost::format("Add residue %1% to pocket.") % (id + 1));
		}
	}
	else
		handle_warning("Given residue ID out of range.");
}

void ProAnalysis::remove_from_pocket(list<size_t> & pocket, size_t id)
{
	if (id < ProE.get_resn())
	{
		if (in_pocket(pocket, id))
		{
			for (list<size_t>::iterator it = pocket.begin(); it != pocket.end(); ++it)
				if (*it == id)
				{
					pocket.erase(it);
					handle_info(boost::format("Remove residue %1% from pocket.") % (id + 1));
				}
		}
		else
			handle_info("Given residue ID not in pocket. Can not erase it.");
	}
	else
		handle_warning("Given residue ID out of range.");
}

void ProAnalysis::gen_pocket(bool has_ligand, list<size_t> &pocket, double cutoff, VectorXd dist2ligand)
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
			handle_warning(boost::format("Given cutoff is too short. Minimum possible cutoff is %1$.2f A.") % dist2ligand.minCoeff());
	}
	else
		handle_warning("Can not find ligand information.");
}

void ProAnalysis::gen_pocket_force(bool & flag, VectorXd &pocket_force, list<size_t> pocket, VectorXd pro_force, VectorXd displacement)
{
	pocket_force = VectorXd::Zero(covariance.rows());

	size_t ndim = pocket.size() * 3;
	MatrixXd X = MatrixXd::Zero(covariance.rows(), ndim);
	VectorXd Y = displacement;
	VectorXd coeff = VectorXd::Zero(ndim);
	//MatrixXd X = MatrixXd::Zero(ndim, ndim);
	//VectorXd Y = VectorXd::Zero(ndim);

	size_t i = 0;
	for (list<size_t>::iterator it = pocket.begin(); it != pocket.end(); ++it)
	{
		/*
		X(3 * i, 3 * i) = covariance(*it * 3, *it * 3);
		X(3 * i, 3 * i + 1) = covariance(*it * 3, *it * 3 + 1);
		X(3 * i, 3 * i + 2) = covariance(*it * 3, *it * 3 + 2);
		X(3 * i + 1, 3 * i) = covariance(*it * 3 + 1, *it * 3);
		X(3 * i + 1, 3 * i + 1) = covariance(*it * 3 + 1, *it * 3 + 1);
		X(3 * i + 1, 3 * i + 2) = covariance(*it * 3 + 1, *it * 3 + 2);
		X(3 * i + 2, 3 * i) = covariance(*it * 3 + 2, *it * 3);
		X(3 * i + 2, 3 * i + 1) = covariance(*it * 3 + 2, *it * 3 + 1);
		X(3 * i + 2, 3 * i + 2) = covariance(*it * 3 + 2, *it * 3 + 2);

		Y(3 * i) = displacement(*it * 3);
		Y(3 * i + 1) = displacement(*it * 3 + 1);
		Y(3 * i + 2) = displacement(*it * 3 + 2);
		*/
		X.col(3 * i) = covariance.col(*it * 3);
		X.col(3 * i + 1) = covariance.col(*it * 3 + 1);
		X.col(3 * i + 2) = covariance.col(*it * 3 + 2);
		
		coeff(i * 3) = pro_force(*it * 3);
		coeff(i * 3 + 1) = pro_force(*it * 3 + 1);
		coeff(i * 3 + 2) = pro_force(*it * 3 + 2);
		
		++i;
	}  
	/*size_t n = coeff.size();
	VectorXd Z = pro_force - coeff;
	std::cout << "verify2£º" << Z << std::endl;
	std::cout << "verify2£º" << n << std::endl;
	*/
 
	switch_LFmethod(coeff, X, Y);

	i = 0;
	for (list<size_t>::iterator it = pocket.begin(); it != pocket.end(); ++it)
	{
		pocket_force(*it * 3) = coeff(i * 3);
		pocket_force(*it * 3 + 1) = coeff(i * 3 + 1);
		pocket_force(*it * 3 + 2) = coeff(i * 3 + 2);
		// add residue will cause disorder ID problem ?
		++i;
	}
	flag = true;
}

void ProAnalysis::gen_pocket_force(bool & flag, VectorXd & pocket_force, VectorXd fixed_force, list<size_t> pocket, list<size_t> fixed_pocket, VectorXd pro_force, VectorXd displacement)
{
	VectorXd equiv_displacement = displacement - covariance * fixed_force / kB / Temp / Navo;
	list<size_t> unfixed_pocket;
	for (list<size_t>::iterator it = pocket.begin(); it != pocket.end(); ++it)
		if (!in_pocket(fixed_pocket, *it))
			unfixed_pocket.push_back(*it);
	gen_pocket_force(flag, pocket_force, unfixed_pocket, pro_force, equiv_displacement);
	pocket_force += fixed_force;
	flag = true;
}

void ProAnalysis::calc_energy_known(bool flag, double &proenergy, double &pocketenergy, double &totenergy, list<size_t> pocket, VectorXd pro_force, MatrixXd distmat, VectorXd displacement)
{
	if (flag)
	{
		VectorXd force = VectorXd::Zero(pro_force.rows());
		for (list<size_t>::iterator it = pocket.begin(); it != pocket.end(); ++it)
		{
			force(*it * 3) = pro_force(*it * 3);
			force(*it * 3 + 1) = pro_force(*it * 3 + 1);
			force(*it * 3 + 2) = pro_force(*it * 3 + 2);
		}
		ArrayXXd distdiffmat = distmat - ProE.get_distmat();
		proenergy = (distdiffmat.pow(2) * ProE.get_kmat()/4).sum();
		pocketenergy = -force.transpose() * displacement;
		//pocketenergy = -force.transpose() * covariance / kB / Temp / Navo * force;
		totenergy = proenergy + pocketenergy;
	}
	else
		handle_warning("Can not calculate energy without pocket force generated. Call \"gen_pocket*_force\" function first.");
}

void ProAnalysis::calc_energy_unknown(bool flag, double &proenergy, double &pocketenergy, double &totenergy, VectorXd pocket_force)
{
	if (flag)
	{
		VectorXd procoord = covariance / kB / Temp / Navo * pocket_force + ProE.get_procoord();
		ArrayXXd distdiffmat = gen_distmat(procoord) - ProE.get_distmat();
		proenergy = (distdiffmat.pow(2) * ProE.get_kmat()/4).sum(); // repeat twice!!
		pocketenergy = -pocket_force.transpose() * covariance / kB / Temp / Navo * pocket_force;
		totenergy = proenergy + pocketenergy;
	}
	else
		handle_warning("Can not calculate energy without pocket force generated. Call \"gen_pocket*_force\" function first.");
}

void ProAnalysis::print_energy_results()
{
	vector<string> buf;
	
	buf.push_back((boost::format("Free energy for binding state structure S: %1$.3f J/mol.") % S_energy).str());
	buf.push_back((boost::format("Pro: %1$.3f J/mol.") % S_proenergy).str());
	buf.push_back((boost::format("Pocket: %1$.3f J/mol.") % S_pocketenergy).str());
	buf.push_back((boost::format("Free energy for allostery state structure A: %1$.3f J/mol.") % A_energy).str());
	buf.push_back((boost::format("Pro: %1$.3f J/mol.") % A_proenergy).str());
	buf.push_back((boost::format("Pocket: %1$.3f J/mol.") % A_pocketenergy).str());
	buf.push_back((boost::format("Free energy for complex state structure AS: %1$.3f J/mol.") % AS_energy).str());
	buf.push_back((boost::format("Pro: %1$.3f J/mol.") % AS_proenergy).str());
	buf.push_back((boost::format("Pocket: %1$.3f J/mol.") % AS_pocketenergy).str());
	buf.push_back((boost::format("Change of free energy: : %1$.3f J/mol.") % ddG).str());
	handle_result("Free energy results: ", buf);
	
	buf.clear();	
	buf.push_back((boost::format("Free energy for binding state structure S: %1$.3f J/mol.") % S_predict_energy).str());
	buf.push_back((boost::format("Pro: %1$.3f J/mol.") % S_predict_proenergy).str());
	buf.push_back((boost::format("Pocket: %1$.3f J/mol.") % S_predict_pocketenergy).str());
	buf.push_back((boost::format("Free energy for allostery state structure A: %1$.3f J/mol.") % A_predict_energy).str());
	buf.push_back((boost::format("Pro: %1$.3f J/mol.") % A_predict_proenergy).str());
	buf.push_back((boost::format("Pocket: %1$.3f J/mol.") % A_predict_pocketenergy).str());
	buf.push_back((boost::format("Free energy for complex state structure AS: %1$.3f J/mol.") % AS_predict_energy).str());
	buf.push_back((boost::format("Pro: %1$.3f J/mol.") % AS_predict_proenergy).str());
	buf.push_back((boost::format("Pocket: %1$.3f J/mol.") % AS_predict_pocketenergy).str());
	buf.push_back((boost::format("Change of free energy: : %1$.3f J/mol.") % ddG_predict).str());
	handle_result("All predict free energy results: ", buf);
}

