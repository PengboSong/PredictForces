#include "ProAnalysis.h"

ProAnalysis::ProAnalysis()
{
	init_container();
}

ProAnalysis::ProAnalysis(Pro apo, Pro binding, Pro allostery, Pro complex)
{
	init_container();

	ProE = apo;
	ProS = binding;
	ProA = allostery;
	ProAS = complex;

	if (!ProE.empty())
	{
		apo_procoord = ProE.get_procoord();
		MinPocket = Minimization(BGDparameters, ProE.get_resn(), apo_procoord, ProE.get_kmat());
		if (!ProS.empty())
		{
			if (ProE.get_resn() != ProS.get_resn())
				Log.handle_message(MSG_ERROR, "Apo state protein and binding state protein do not have equal residue numbers.");

			preprocess(POCKETS);
		}
		if (!ProA.empty())
		{
			if (ProE.get_resn() != ProA.get_resn())
				Log.handle_message(MSG_ERROR, "Apo state protein and allostery state protein do not have equal residue numbers.");

			preprocess(POCKETA);
		}
		if (!ProAS.empty())
		{
			if (ProE.get_resn() != ProAS.get_resn())
				Log.handle_message(MSG_ERROR, "Apo state protein and complex state protein do not have equal residue numbers.");

			preprocess(POCKETAS);
		}
	}
}

ProAnalysis::~ProAnalysis()
{
}

void ProAnalysis::init_container()
{
	PocketInfo pocketS, pocketA, pocketAS;
	pockets[POCKETS] = pocketS;
	pockets[POCKETA] = pocketA;
	pockets[POCKETAS] = pocketAS;

	ProInfo proinfoS, proinfoA, proinfoAS;
	proinfos[POCKETS] = proinfoS;
	proinfos[POCKETA] = proinfoA;
	proinfos[POCKETAS] = proinfoAS;
}

void ProAnalysis::preprocess(Pockets m)
{
	Pro* targetpro = &ProE;
	std::string keyword;
	switch (m)
	{
	case POCKETS:
		*targetpro = ProS;
		keyword = "binding";
		break;
	case POCKETA:
		*targetpro = ProA;
		keyword = "allostery";
		break;
	case POCKETAS:
		*targetpro = ProAS;
		keyword = "complex";
		break;
	}

	/* <--- Protein ---> */

	pro(m).empty = false;
	pro(m).withligand = targetpro->has_ligand();
	pro(m).dist2ligand = targetpro->get_dist2ligand();
	pro(m).coord = targetpro->get_procoord();

	pro(m).fitcoord = fitting(apo_procoord, pro(m).coord);
	Log.handle_message(MSG_INFO, "Fitting process succeed.");
	pro(m).displacement = calc_displacement(apo_procoord, pro(m).coord);
	Log.handle_message(MSG_INFO, "Calculating displacement succeed.");
	pro(m).rmsd = calc_rmsd(pro(m).displacement);
	Log.handle_message(MSG_RESULT, boost::format("RMSD between apo state and %1% state from PDB file: %2$.4f A.") % keyword % pro(m).rmsd);
	
	gen_pocket(m, pocket_cutoff);
	show_pocket(m);

	/* <--- Pocket ---> */

	size_t pocketn = pocket_members(m).size() * 3;
	VectorXd fixed_apo = VectorXd::Zero(pocketn), fixed_holo = VectorXd::Zero(pocketn);

	grep_pocket_coord(fixed_apo, apo_procoord, pocket_members(m));
	grep_pocket_coord(fixed_holo, pro(m).coord, pocket_members(m));

	pocket(m).coord = fitting(fixed_apo, fixed_holo);
	Log.handle_message(MSG_INFO, "Fitting pocket process succeed.");
	pocket(m).displacement = calc_displacement(fixed_apo, pocket(m).coord);
	Log.handle_message(MSG_INFO, "Calculating pocket displacement succeed.");
	pocket(m).rmsd = calc_rmsd(pocket(m).displacement);
	Log.handle_message(MSG_RESULT, boost::format("RMSD from %1% state PDB file: %2$.4f A.") % keyword % pocket(m).rmsd);

	pro(m).preprocess = true;
}

void ProAnalysis::interactive_pocket(Pockets m)
{
	std::string buf, cmd, label, mode;
	std::vector<std::string> para;

	switch (m)
	{
	case POCKETS:
		mode = "S";
		break;
	case POCKETA:
		mode = "A";
		break;
	case POCKETAS:
		mode = "AS";
		break;
	default:
		mode = "Err";
	}

	while (true)
	{
		std::cout << '[' << mode << ']' << " >>> ";
		std::getline(std::cin, buf);
		boost::algorithm::trim(buf);
		split(para, buf, boost::is_any_of(" "));
		cmd = para[0];
		if (cmd == "back")
			break;
		else if (cmd == "add")
		{
			size_t resid = 0;
			bool convert_flag = true;
			try
			{
				resid = boost::lexical_cast<size_t>(para[1]) - 1;
			}
			catch (boost::bad_lexical_cast)
			{
				convert_flag = false;
				HandleMessage::print(MSG_EMPTY, "No residue ID given. Please enter again.");
			}

			if (convert_flag)
				add_to_pocket(m, resid);
		}
		else if (cmd == "del")
		{
			size_t resid = 0;
			bool convert_flag = true;
			try
			{
				resid = boost::lexical_cast<size_t>(para[1]) - 1;
			}
			catch (boost::bad_lexical_cast)
			{
				convert_flag = false;
				HandleMessage::print(MSG_EMPTY, "No residue ID given. Please enter again.");
			}

			if (convert_flag)
				remove_from_pocket(m, resid);
		}
		else if (cmd == "gen-pocket")
		{
			double cutoff = 0.0;
			bool convert_flag = true;
			try
			{
				cutoff = boost::lexical_cast<double>(para[1]);
			}
			catch (boost::bad_lexical_cast)
			{
				convert_flag = false;
				HandleMessage::print(MSG_EMPTY, "Illegal cutoff length given. Please enter again.");
			}

			if (convert_flag)
				gen_pocket(m, cutoff);
		}
		else if (cmd == "show")
		{
			show_pocket(m);
		}
		else if (cmd == "test")
		{
			test_pocket(m);
		}
		else if (cmd == "gen-force")
		{
			gen_pocket_force(m);
		}
		else if (cmd == "show-force")
		{
			show_pocket_force(m);
		}
		else if (cmd == "origin-force")
		{
			show_pro_pocket_force(m);
		}
		else if (cmd == "all-origin-force")
		{
			show_pro_all_force(m);
		}
		else
		{
			HandleMessage::print(MSG_EMPTY, "Unknown command.");
		}
		cmd.clear();
	}
}

void ProAnalysis::interactive()
{
	std::string buf, cmd;
	std::vector<std::string> para;
	while (true)
	{
		std::cout << ">>> ";
		std::getline(std::cin, buf);
		boost::algorithm::trim(buf);
		split(para, buf, boost::is_any_of(" "));
		cmd = para[0];
		if (cmd == "exit")
			break;
		else if (cmd == "pocketS")
			interactive_pocket(POCKETS);
		else if (cmd == "pocketA")
			interactive_pocket(POCKETA);
		else if (cmd == "pocketAS")
			interactive_pocket(POCKETAS);
		else if (cmd == "energy")
			gen_free_energy();
		else if (cmd == "LFmethod")
			show_LFmethod_detail();
		else
			HandleMessage::print(MSG_EMPTY, "Unknown command.");
		cmd.clear();
	}
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

void ProAnalysis::show_LFmethod_detail()
{
	switch (LFMETHOD)
	{
	case NormalEquation:
		// Normal Equation
		Log.handle_message(MSG_INFO, "Multiple Linear Fitting using normal equation.");
		break;
	case BatchGradientDescent:
		// Batch Gradient Descent
		Log.handle_message(MSG_INFO, "Multiple Linear Fitting using batch gradient descent.");
		Log.handle_message(MSG_INFO, boost::format("Using learning step %1$.2e") % BGDparameters.learning_rate);
		Log.handle_message(MSG_INFO, boost::format("Using convergence %1$.2e") % BGDparameters.convergence);
		Log.handle_message(MSG_INFO, boost::format("Using maximum iteration times %1%") % BGDparameters.niteration);
		break;
	default:
		Log.handle_message(MSG_WARNING, "Invalid LF method mode.");
	}
}

void ProAnalysis::gen_pocket_force(Pockets m)
{
	if (pro(m).preprocess)
	{
		if (m == POCKETAS)
		{
			if (pocket(POCKETS).access_force)
			{
				gen_pocket_force(
					pocket(POCKETAS).access_force,
					pocket(POCKETAS).force,
					pocket(POCKETS).force,
					pocket(POCKETAS).members,
					pocket(POCKETS).members,
					pro(POCKETAS).force,
					pro(POCKETAS).displacement
				);
			}
			else if (pocket(POCKETA).access_force)
			{
				gen_pocket_force(
					pocket(POCKETAS).access_force,
					pocket(POCKETAS).force,
					pocket(POCKETA).force,
					pocket(POCKETAS).members,
					pocket(POCKETA).members,
					pro(POCKETAS).force,
					pro(POCKETAS).displacement
				);
			}
			else
			{
				gen_pocket_force(
					pocket(POCKETAS).access_force,
					pocket(POCKETAS).force,
					pocket(POCKETAS).members,
					pro(POCKETAS).force,
					pro(POCKETAS).displacement
				);
			}
		}
		else
		{
			gen_pocket_force(
				pocket(m).access_force,
				pocket(m).force,
				pocket(m).members,
				pro(m).force,
				pro(m).displacement
			);
		}
	}
}

void ProAnalysis::gen_free_energy()
{
	if (pro(POCKETS).preprocess && pro(POCKETA).preprocess)
	{
		energy_minimization(FixPocket, POCKETS);
		check_EM_diff(POCKETS);

		energy_minimization(FixPocket, POCKETA);
		check_EM_diff(POCKETA);

		pocket(POCKETAS).force = pocket(POCKETS).force + pocket(POCKETA).force;
		energy_minimization(OnlyForce, POCKETS);

        ddG = pro(POCKETAS).G.total - pro(POCKETS).G.total - pro(POCKETA).G.total;

		print_energy_results();
	}
	else if ((pro(POCKETS).preprocess || pro(POCKETA).preprocess) && pro(POCKETAS).preprocess)
	{
		Pockets s;	// POCKET identifier for single pocket
		Pockets as; // POCKET identifier for another pocket
		if (pro(POCKETS).preprocess)
		{
			s = POCKETS;
			as = POCKETA;
		}
		else
		{
			s = POCKETA;
			as = POCKETS;
		}

		energy_minimization(FixPocket, s);
		check_EM_diff(s);
		pocket(POCKETAS).force = pocket(s).force;

		energy_minimization(FixPocketWithForce, POCKETAS);
		check_EM_diff(POCKETAS);

		size_t pocketn = pocket_members(s).size() * 3;
		VectorXd pocketcoord = VectorXd::Zero(pocketn), pocketcoord_apo = VectorXd::Zero(pocketn);
		grep_pocket_coord(pocketcoord, pro(POCKETAS).equilibrium_coord, pocket_members(s));
		grep_pocket_coord(pocketcoord_apo, apo_procoord, pocket_members(s));

		double rmsd_pocket = calc_rmsd(pocketcoord, fitting(pocketcoord_apo, pocketcoord));
		Log.handle_message(MSG_RESULT, boost::format("RMSD of pocket coordinates with domain coordinates in apo state structure: %1$.4f A.") % rmsd_pocket);

		pocket(as).force = pocket(POCKETAS).force - pocket(s).force;
		energy_minimization(OnlyForce, as);
		
		ddG = pro(POCKETAS).G.total - pro(s).G.total - pro(as).G.total;
		print_energy_results();
	}
	else
		Log.handle_message(MSG_ERROR, "Lack necessary protein information.");
}

void ProAnalysis::write_matrix(MatrixXd mat, std::string writepath)
{
	std::ofstream matf(writepath, std::ios::out);
	if (matf.is_open())
	{
		matf << mat.format(CleanFmt);
		matf.close();
		Log.handle_message(MSG_INFO, boost::format("Matrix has been written to %1%.") % writepath);
	}
	else
		Log.handle_message(MSG_WARNING, boost::format("Can not open file %1%.") % writepath);
}

void ProAnalysis::write_matrix_binary(MatrixXd mat, std::string writepath)
{
	std::ofstream matf(writepath, std::ios::out | std::ios::binary);
	if (matf.is_open())
	{
		size_t Msize = mat.size();
		size_t Mrown = mat.rows();
		size_t Mcoln = mat.cols();
		double* M = new double[Msize];
		for (size_t i = 0; i < Mrown; ++i)
			for (size_t j = 0; j < Mcoln; ++j)
				M[i * Mcoln + j] = mat(i, j);
		matf.write((char *)&Msize, sizeof(size_t));
		matf.write((char *)&Mrown, sizeof(size_t));
		matf.write((char *)&Mcoln, sizeof(size_t));
		matf.write((char *)&M[0], Msize * sizeof(double));
		matf.close();
		delete[] M;
		Log.handle_message(MSG_INFO, boost::format("Matrix (binary data) has been written to %1%.") % writepath);
	}
	else
		Log.handle_message(MSG_WARNING, boost::format("Can not open file %1%.") % writepath);
}

void ProAnalysis::read_matrix_binary(MatrixXd & mat, std::string fpath)
{
	std::ifstream matf(fpath, std::ios::in | std::ios::binary);
	if (matf.is_open())
	{
		size_t Msize = 0;
		size_t Mrown = 0;
		size_t Mcoln = 0;
		matf.read((char *)&Msize, sizeof(size_t));
		matf.read((char *)&Mrown, sizeof(size_t));
		matf.read((char *)&Mcoln, sizeof(size_t));
		double* M = new double[Msize];
		matf.read((char *)&M[0], Msize * sizeof(double));
		mat = MatrixXd::Zero(Mrown, Mcoln);
		for (size_t i = 0; i < Mrown; ++i)
			for (size_t j = 0; j < Mcoln; ++j)
				mat(i, j) = M[i * Mcoln + j];
		matf.close();
		delete[] M;
	}
	else
		Log.handle_message(MSG_WARNING, boost::format("Can not open file %1%.") % fpath);
}

std::string ProAnalysis::write_path(std::string fname)
{
	filesys::path workdir = workdirPath;
	return (workdir / fname).string();
}

void ProAnalysis::switch_LFmethod(VectorXd & coeff, MatrixXd X, VectorXd Y)
{
	switch (LFMETHOD)
	{
	case NormalEquation:
		normal_equation(coeff, X, Y);
		break;
	case BatchGradientDescent:
		BGD(coeff, X, Y, BGDparameters);
		break;
	default:
		BGD(coeff, X, Y, BGDparameters);
	}
}

void ProAnalysis::energy_minimization(EMType em, Pockets m)
{
	switch (em)
	{
	case FixPocket:
		MinPocket.single_pocket(pocket(m).members, pocket(m).coord);
		break;
	case FixPocketWithForce:
		MinPocket.single_pocket(pocket(m).members, pocket(m).force, pocket(m).coord);
		break;
	case OnlyForce:
		MinPocket.single_pocket(pocket(m).force);
		break;
	default:
		Log.handle_message(MSG_WARNING, "Energy minimization task type is not specified.");
	}
	if (!MinPocket.complete())
		return;
	pocket(m).access_force = MinPocket.convergence();
	pocket(m).force = MinPocket.pocket_force();
	pro(m).equilibrium_coord = MinPocket.equilibrium_coord();
	pro(m).G = MinPocket.energy();
}

double ProAnalysis::calc_model_rmsd(bool flag, VectorXd pocket_force, VectorXd refcoord)
{
	if (flag)
	{
		mprocoord = (covariance / kB / Temp / Navo) * pocket_force  + ProE.get_procoord();
		VectorXd fitmprocoord = fitting(refcoord, mprocoord);
		return calc_rmsd(refcoord, fitmprocoord);
	}
	else
	{
		Log.handle_message(MSG_WARNING, "Can not calculate model RMSD without pocket force generated. Call \"gen_pocket*_force\" function first.");
		return 0.0;
	}
}

double ProAnalysis::calc_model_correlation(bool flag, VectorXd pocket_force, VectorXd displacement)
{
	if (flag)
	{
		mprocoord = (covariance / kB / Temp / Navo) * pocket_force + ProE.get_procoord();
		VectorXd new_displacement = calc_displacement(ProE.get_procoord(), mprocoord);

		if (displacement.size() == new_displacement.size())
		{
			double xBar = calc_average(new_displacement);
			double yBar = calc_average(displacement);
			double varX = 0;
			double varY = 0;
			double SSR = 0;
			double SST = 0;

			size_t vlen = displacement.size() / 3;

			VectorXd X = VectorXd::Zero(vlen);
			VectorXd Y = VectorXd::Zero(vlen);
			VectorXd diffX = VectorXd::Zero(vlen);
			VectorXd diffY = VectorXd::Zero(vlen);

			for (size_t i = 0; i < vlen; ++i)
			{
				X(i) = sqrt(pow(new_displacement(i * 3), 2) + pow(new_displacement(i * 3 + 1), 2) + pow(new_displacement(i * 3 + 2), 2));
				Y(i) = sqrt(pow(displacement(i * 3), 2) + pow(displacement(i * 3 + 1), 2) + pow(displacement(i * 3 + 2), 2));
				diffX(i) = X(i) - xBar;
				diffY(i) = Y(i) - yBar;
				SSR += (diffX(i) * diffY(i));
				varX += pow(diffX(i), 2);
				varY += pow(diffY(i), 2);
			}
				
			SST = sqrt(varX * varY);
			return SSR / SST;
		}
		else
			return 0.0;
	}
	else
	{
		Log.handle_message(MSG_WARNING, "Can not calculate model correlation without pocket force generated. Call \"gen_pocket*_force\" function first.");
		return 0.0;
	}
}

double ProAnalysis::calc_correlation(VectorXd displacement, VectorXd new_displacement)
{
	if (displacement.size() == new_displacement.size())
	{
		double xBar = calc_average(new_displacement);
		double yBar = calc_average(displacement);
		double varX = 0;
		double varY = 0;
		double SSR = 0;
		double SST = 0;

		size_t vlen = displacement.size() / 3;

		VectorXd X = VectorXd::Zero(vlen);
		VectorXd Y = VectorXd::Zero(vlen);
		VectorXd diffX = VectorXd::Zero(vlen);
		VectorXd diffY = VectorXd::Zero(vlen);

		for (size_t i = 0; i < vlen; ++i)
		{
			X(i) = sqrt(pow(new_displacement(i * 3), 2) + pow(new_displacement(i * 3 + 1), 2) + pow(new_displacement(i * 3 + 2), 2));
			Y(i) = sqrt(pow(displacement(i * 3), 2) + pow(displacement(i * 3 + 1), 2) + pow(displacement(i * 3 + 2), 2));
			diffX(i) = X(i) - xBar;
			diffY(i) = Y(i) - yBar;
			SSR += (diffX(i) * diffY(i));
			varX += pow(diffX(i), 2);
			varY += pow(diffY(i), 2);
		}

		SST = sqrt(varX * varY);
		return SSR / SST;
	}
	else
		return 0.0;
}



void ProAnalysis::show_pocket(PocketList pocket)
{
	std::string buf = "Pocket residues: ";
	for (PocketList::iterator it = pocket.begin(); it != pocket.end(); ++it)
	{
		buf += std::to_string(*it + 1) + " ";
	}

	Log.handle_message(MSG_RESULT, buf);
}
void ProAnalysis::test_pocket(bool flag, bool info, PocketList pocket, VectorXd displacement, VectorXd pocket_force, VectorXd refcoord)
{
	if (!pocket.empty() && flag)
	{
		show_pocket(pocket);
		if (info)
		{
			Log.handle_message(
				MSG_RESULT,
				boost::format("RMSD between real structure and structure calculated according to current pocket: %1% A.") % calc_model_rmsd(flag, pocket_force, refcoord)
			);
			Log.handle_message(
				MSG_RESULT,
				boost::format("Pearson between real structure and structure calculated according to current pocket: %1% A.") % calc_model_correlation(flag, displacement, pocket_force)
			);
		}
		else
			Log.handle_message(MSG_WARNING, "Lack necessary protein information.");
	}
	else
		Log.handle_message(MSG_WARNING, "The binding pocket domain is not specificed.");
}

void ProAnalysis::show_pocket_force(bool access, PocketList pocket_members, VectorXd pocket_force)
{
	if (access)
	{
		Log.handle_message(MSG_RESULT, "POCKET FORCE (PREDICT)");
		Log.handle_message(MSG_EMPTY, "*---*---*---*---*---*");
		for (PocketList::const_iterator it = pocket_members.cbegin(); it != pocket_members.cend(); ++it)
		{
			Vector3d resforce = Vector3d::Zero();
			resforce << pocket_force(*it * 3), pocket_force(*it * 3 + 1), pocket_force(*it * 3 + 2);
			Log.handle_message(MSG_EMPTY, boost::format("RES %1% FORCE %2$.4f") % (*it + 1) % calc_norm(resforce));
		}
		Log.handle_message(MSG_EMPTY, "*---*---*---*---*---*");
	}
	else
		Log.handle_message(MSG_WARNING, "The binding pocket domain is not specificed.");
}

void ProAnalysis::show_pro_pocket_force(PocketList pocket, VectorXd pro_force)
{
	Log.handle_message(MSG_RESULT, "POCKET FORCE (STRUCTURE)");
	Log.handle_message(MSG_EMPTY, "*---*---*---*---*---*");
	for (PocketList::iterator it = pocket.begin(); it != pocket.end(); ++it)
	{
		Vector3d resforce = Vector3d::Zero();
		resforce << pro_force(*it * 3), pro_force(*it * 3 + 1), pro_force(*it * 3 + 2);
		Log.handle_message(MSG_EMPTY, boost::format("RES %1% FORCE %2$.4f") % (*it + 1) % calc_norm(resforce));
	}
	Log.handle_message(MSG_EMPTY, "*---*---*---*---*---*");
}

void ProAnalysis::show_pro_all_force(VectorXd pro_force)
{
	std::vector<std::string> buf;
	std::vector<std::pair<size_t, double>> forces;
	for (size_t i = 0; i < ProE.get_resn(); ++i)
	{
		Vector3d resforce = Vector3d::Zero();
		resforce << pro_force(i * 3), pro_force(i * 3 + 1), pro_force(i * 3 + 2);
		forces.push_back(std::make_pair(i, calc_norm(resforce)));
	}
	std::sort(forces.begin(), forces.end(), [](std::pair<size_t, double> a, std::pair<size_t, double> b) {
		return a.second > b.second;
	});

	Log.handle_message(MSG_RESULT, "PROTEIN FORCE");
	Log.handle_message(MSG_EMPTY, "*---*---*---*---*---*");
	for (size_t i = 0; i < ProE.get_resn(); ++i)
		Log.handle_message(MSG_EMPTY, boost::format("RES %1% FORCE %2$.4f") % (forces[i].first + 1) % forces[i].second);
	Log.handle_message(MSG_EMPTY, "*---*---*---*---*---*");
}

bool ProAnalysis::in_pocket(PocketList pocket, size_t id)
{
	bool find_element_flag = false;
	for (PocketList::iterator it = pocket.begin(); it != pocket.end(); ++it)
		if (*it == id)
			find_element_flag = true;
	return find_element_flag;
}

void ProAnalysis::add_to_pocket(PocketList & pocket, size_t id)
{
	if (id < ProE.get_resn())
	{
		if (in_pocket(pocket, id))
			Log.handle_message(MSG_INFO, "Given residue ID already in pocket.");
		else
		{
			pocket.push_back(id);
			Log.handle_message(MSG_INFO, boost::format("Add residue %1% to pocket.") % (id + 1));
		}
	}
	else
		Log.handle_message(MSG_WARNING, "Given residue ID out of range.");
}

void ProAnalysis::remove_from_pocket(PocketList & pocket, size_t id)
{
	if (id < ProE.get_resn())
	{
		if (in_pocket(pocket, id))
		{
			for (PocketList::iterator it = pocket.begin(); it != pocket.end(); ++it)
				if (*it == id)
				{
					pocket.erase(it);
					Log.handle_message(MSG_INFO, boost::format("Remove residue %1% from pocket.") % (id + 1));
				}
		}
		else
			Log.handle_message(MSG_INFO, "Given residue ID not in pocket. Can not erase it.");
	}
	else
		Log.handle_message(MSG_WARNING, "Given residue ID out of range.");
}

void ProAnalysis::gen_pocket(bool has_ligand, PocketList & pocket, double cutoff, VectorXd dist2ligand)
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
			Log.handle_message(
				MSG_WARNING,
				boost::format("Given cutoff is too short. Minimum possible cutoff is %1$.2f A.") % dist2ligand.minCoeff()
			);
	}
	else
		Log.handle_message(MSG_WARNING, "Can not find ligand information.");
}

void ProAnalysis::gen_pocket_force(bool & flag, VectorXd & pocket_force, PocketList pocket, VectorXd pro_force, VectorXd displacement)
{
	pocket_force = VectorXd::Zero(covariance.rows());  
	size_t ndim = pocket.size() * 3;
	
	MatrixXd X = MatrixXd::Zero(covariance.rows(), ndim);  // consider all residues for cost function.
	VectorXd Y = displacement;
	VectorXd coeff = VectorXd::Zero(ndim);
	grep_pocket_coord(coeff, pro_force, pocket);

	size_t i = 0;
	for (PocketList::iterator it = pocket.begin(); it != pocket.end(); ++it)  
	{
		X.col(3 * i) = covariance.col(*it * 3);
		X.col(3 * i + 1) = covariance.col(*it * 3 + 1);
		X.col(3 * i + 2) = covariance.col(*it * 3 + 2);
		++i;
	}

	X /= (Navo * Temp * kB);

	switch_LFmethod(coeff, X, Y);

	modify_pocket_coord(pocket_force, coeff, pocket);

	flag = true;
}

void ProAnalysis::gen_pocket_force(bool & flag, VectorXd & pocket_force, VectorXd fixed_force, PocketList pocket, PocketList fixed_pocket, VectorXd pro_force, VectorXd displacement)
{
	VectorXd equiv_displacement = displacement - covariance * fixed_force / kB / Temp / Navo;
	PocketList unfixed_pocket;
	for (PocketList::iterator it = pocket.begin(); it != pocket.end(); ++it)
		if (!in_pocket(fixed_pocket, *it))
			unfixed_pocket.push_back(*it);
	gen_pocket_force(flag, pocket_force, unfixed_pocket, pro_force, equiv_displacement);
	pocket_force += fixed_force;
	flag = true;
}

void ProAnalysis::check_EM_diff(Pockets m)
{
	double rmsd_em = calc_rmsd(pro(m).equilibrium_coord, fitting(pro(m).equilibrium_coord, pro(m).coord));
	Log.handle_message(
		MSG_RESULT,
		boost::format("RMSD of coordinates between PDB file structure and minimized structure: %1$.4f A.") % rmsd_em
	);

	VectorXd displacement = fitting(apo_procoord, pro(m).coord) - apo_procoord;
	VectorXd displacement_em = pro(m).equilibrium_coord - apo_procoord;
	double pearson = calc_correlation(displacement, displacement_em);
	Log.handle_message(
		MSG_RESULT,
		boost::format("Pearson between PDB file structure and minimized structure: %1$.4f") % pearson
	);
}

void ProAnalysis::print_energy_results()
{
	Log.handle_message(MSG_RESULT, "Energy results:");
	Log.handle_message(MSG_EMPTY, "*---*---*---*---*---*");
	Log.handle_message(
		MSG_EMPTY,
		boost::format("Free energy for binding state structure S: %1$.3f Kcal/mol.") % pro(POCKETS).G.total
	);
	Log.handle_message(
		MSG_EMPTY,
		boost::format("Pro: %1$.3f Kcal/mol.") % pro(POCKETS).G.pro
	);
	Log.handle_message(
		MSG_EMPTY,
		boost::format("Pocket: %1$.3f Kcal/mol.") % pro(POCKETS).G.pocket
	);
	Log.handle_message(
		MSG_EMPTY,
		boost::format("Free energy for allostery state structure A: %1$.3f Kcal/mol.") % pro(POCKETA).G.total
	);
	Log.handle_message(
		MSG_EMPTY,
		boost::format("Pro: %1$.3f Kcal/mol.") % pro(POCKETA).G.pro
	);
	Log.handle_message(
		MSG_EMPTY,
		boost::format("Pocket: %1$.3f Kcal/mol.") % pro(POCKETA).G.pocket
	);
	Log.handle_message(
		MSG_EMPTY,
		boost::format("Free energy for complex state structure AS: %1$.3f Kcal/mol.") % pro(POCKETAS).G.total
	);
	Log.handle_message(
		MSG_EMPTY,
		boost::format("Pro: %1$.3f Kcal/mol.") % pro(POCKETAS).G.pro
	);
	Log.handle_message(
		MSG_EMPTY,
		boost::format("Pocket: %1$.3f Kcal/mol.") % pro(POCKETAS).G.pocket
	);
	Log.handle_message(
		MSG_EMPTY,
		boost::format("Change of free energy: %1$.3f Kcal/mol.") % ddG
	);
	Log.handle_message(MSG_EMPTY, "*---*---*---*---*---*");
}
