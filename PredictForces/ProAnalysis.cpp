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
		if (!ProS.empty())
		{
			if (ProE.get_resn() != ProS.get_resn())
				handle_message(MSG_ERROR, "Apo state protein and binding state protein do not have equal residue numbers.");

			preprocess(POCKETS);
		}
		if (!ProA.empty())
		{
			if (ProE.get_resn() != ProA.get_resn())
				handle_message(MSG_ERROR, "Apo state protein and allostery state protein do not have equal residue numbers.");

			preprocess(POCKETA);
		}
		if (!ProAS.empty())
		{
			if (ProE.get_resn() != ProAS.get_resn())
				handle_message(MSG_ERROR, "Apo state protein and complex state protein do not have equal residue numbers.");

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

	ApoProInfo apo_proinfoS, apo_proinfoA, apo_proinfoAS;
	apo_proinfos[POCKETS] = apo_proinfoS;
	apo_proinfos[POCKETA] = apo_proinfoA;
	apo_proinfos[POCKETAS] = apo_proinfoAS;
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

	pro(m).empty = false;
	pro(m).withligand = targetpro->has_ligand();
	pro(m).dist2ligand = targetpro->get_dist2ligand();
	pro(m).procoord = targetpro->get_procoord();
	gen_pocket(m, pocket_cutoff);
	show_pocket(m);

	size_t ndim = pocket_members(m).size() * 3;
	VectorXd fixed_apo = VectorXd::Zero(ndim);
	VectorXd fixed_holo = VectorXd::Zero(ndim);
	VectorXd init_coord = ProE.get_procoord();
	VectorXd holo_coord = targetpro->get_procoord();

	size_t i = 0;
	for (PocketList::iterator it = pocket_members(m).begin(); it != pocket_members(m).end(); ++it)
	{
		fixed_apo(i * 3) = init_coord(*it * 3);
		fixed_apo(i * 3 + 1) = init_coord(*it * 3 + 1);
		fixed_apo(i * 3 + 2) = init_coord(*it * 3 + 2);

		fixed_holo(i * 3) = holo_coord(*it * 3);
		fixed_holo(i * 3 + 1) = holo_coord(*it * 3 + 1);
		fixed_holo(i * 3 + 2) = holo_coord(*it * 3 + 2);

		++i;
	}

	pro(m).fitprocoord = fitting(fixed_apo, fixed_holo);
	handle_message(MSG_INFO, "Fitting process succeed.");

	apo_pro(m).displacement = calc_displacement(fixed_apo, pro(m).fitprocoord);
	handle_message(MSG_INFO, "Calculating displacement succeed.");

	apo_pro(m).rmsd = calc_rmsd(apo_pro(m).displacement);
	handle_message(MSG_RESULT, boost::format("RMSD from %1% state PDB file: %2$.4f A.") % keyword % apo_pro(m).rmsd);

	apo_pro(m).force = hessian * apo_pro(m).displacement;
	handle_message(MSG_INFO, "Calculating force succeed.");

	apo_pro(m).distdiff = gen_differ(targetpro->get_distmat(), ProE.get_distmat());
	handle_message(MSG_INFO, "Constructing distance difference matrix succeed.");

	apo_pro(m).preprocess = true;
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
				handle_message(MSG_EMPTY, "No residue ID given. Please enter again.");
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
				handle_message(MSG_EMPTY, "No residue ID given. Please enter again.");
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
				handle_message(MSG_EMPTY, "Illegal cutoff length given. Please enter again.");
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
			handle_message(MSG_EMPTY, "Unknown command.");
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
			handle_message(MSG_EMPTY, "Unknown command.");
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
		handle_message(MSG_INFO, "Multiple Linear Fitting using normal equation.");
		break;
	case BatchGradientDescent:
		// Batch Gradient Descent
		handle_message(MSG_INFO, "Multiple Linear Fitting using batch gradient descent.");
		handle_message(MSG_INFO, boost::format("Using learning step %1$.2e") % BGDparameters.learning_rate);
		handle_message(MSG_INFO, boost::format("Using convergence %1$.2e") % BGDparameters.convergence);
		handle_message(MSG_INFO, boost::format("Using maximum iteration times %1%") % BGDparameters.niteration);
		break;
	default:
		handle_message(MSG_WARNING, "Invalid LF method mode.");
	}
}

void ProAnalysis::gen_pocket_force(Pockets m)
{
	if (apo_pro(m).preprocess)
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
					apo_pro(POCKETAS).force,
					apo_pro(POCKETAS).displacement
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
					apo_pro(POCKETAS).force,
					apo_pro(POCKETAS).displacement
				);
			}
			else
			{
				gen_pocket_force(
					pocket(POCKETAS).access_force,
					pocket(POCKETAS).force,
					pocket(POCKETAS).members,
					apo_pro(POCKETAS).force,
					apo_pro(POCKETAS).displacement
				);
			}
		}
		else
		{
			gen_pocket_force(
				pocket(m).access_force,
				pocket(m).force,
				pocket(m).members,
				apo_pro(m).force,
				apo_pro(m).displacement
			);
		}
	}
}

void ProAnalysis::minimization_calc_energy(bool & flag, Pockets m)
{
	minimization(
		flag,
		pocket(m).members,
		pro(m).fitprocoord,
		apo_pro(m).equilibrium_coord,
		pocket(m).force,
		BGDparameters
	);
	calc_energy_unknown(
		pocket(m).access_force,
		pro(m).G,
		pocket(m).force,
		apo_pro(m).equilibrium_coord
	);
}

void ProAnalysis::equilibrium_coord_rmsd(Pockets m)
{

	VectorXd displacement_complex = apo_pro(m).equilibrium_coord - fitting(apo_pro(m).equilibrium_coord, pro(m).procoord);
	double rmsd_complex = calc_rmsd(displacement_complex);
	handle_message(MSG_RESULT, boost::format("RMSD with complex state PDB file: %1$.4f A.") % rmsd_complex);

	VectorXd displacement_apo = fitting(apo_procoord, apo_pro(m).equilibrium_coord) - apo_procoord;
	double rmsd_apo = calc_rmsd(displacement_apo);
	handle_message(MSG_RESULT, boost::format("RMSD with apo state PDB file: %1$.4f A.") % rmsd_apo);
}

void ProAnalysis::gen_free_energy()
{
	if (apo_pro(POCKETS).preprocess && apo_pro(POCKETA).preprocess)
	{
		bool minimization_status = false;
		
		minimization_calc_energy(minimization_status, POCKETS);

		equilibrium_coord_rmsd(POCKETS);

		minimization_calc_energy(minimization_status, POCKETA);

		equilibrium_coord_rmsd(POCKETA);

		pocket(POCKETAS).force = pocket(POCKETS).force + pocket(POCKETA).force;
		
		minimization(minimization_status, apo_pro(POCKETAS).equilibrium_coord, pocket(POCKETAS).force, BGDparameters);
		calc_energy_unknown(pocket(POCKETAS).access_force, pro(POCKETAS).G, pocket(POCKETAS).force, apo_pro(POCKETAS).equilibrium_coord);

        ddG = pro(POCKETAS).G.total - pro(POCKETS).G.total - pro(POCKETA).G.total;

		print_energy_results();
	}
	else if ((apo_pro(POCKETS).preprocess || apo_pro(POCKETA).preprocess) && apo_pro(POCKETAS).preprocess)
	{
		Pockets s;	// POCKET identifier for single pocket
		Pockets as; // POCKET identifier for another pocket
		if (apo_pro(POCKETS).preprocess)
		{
			s = POCKETS;
			as = POCKETA;
		}
		else
		{
			s = POCKETA;
			as = POCKETS;
		}

		bool minimization_status = false;

		minimization_calc_energy(minimization_status, s);

		minimization_calc_energy(minimization_status, POCKETAS);
		
		size_t i = 0, len = pocket_members(s).size() * 3;
		VectorXd new_pocket_coord = VectorXd::Zero(len);
		VectorXd new_apo_pocket_coord = VectorXd::Zero(len);
		grep_pocket_coord(new_pocket_coord, apo_pro(POCKETAS).equilibrium_coord, pocket_members(s));
		grep_pocket_coord(new_apo_pocket_coord, apo_procoord, pocket_members(s));

		double gall = calc_rmsd(new_pocket_coord, fitting(new_apo_pocket_coord, new_pocket_coord));

		pocket(as).force = pocket(POCKETAS).force - pocket(s).force;

		minimization(minimization_status, apo_pro(as).equilibrium_coord, pocket(as).force, BGDparameters);
		calc_energy_unknown(pocket(as).access_force, pro(as).G, pocket(as).force, apo_pro(as).equilibrium_coord);

		ddG = pro(POCKETAS).G.total - pro(s).G.total - pro(as).G.total;

		VectorXd s_pocket_displacement = apo_pro(POCKETS).equilibrium_coord - fitting(apo_pro(s).equilibrium_coord, pro(s).procoord);
		double s_pocket_rmsd = calc_rmsd(s_pocket_displacement);
		handle_message(MSG_RESULT, boost::format("RMSD from complex state PDB file: %1$.4f A.") % s_pocket_rmsd);

		VectorXd AS_displacement = apo_pro(POCKETS).equilibrium_coord - fitting(apo_pro(POCKETAS).equilibrium_coord, pro(POCKETAS).procoord);
		double AS_rmsd = calc_rmsd(AS_displacement);
		handle_message(MSG_RESULT, boost::format("RMSD from complex state PDB file: %1$.4f A.") % AS_rmsd);

		double s_pocket_pearson = calc_correlation(
			fitting(apo_procoord, pro(s).procoord) - apo_procoord,
			apo_pro(POCKETAS).equilibrium_coord - apo_procoord
		);
		handle_message(
			MSG_RESULT,
			boost::format("Pearson between real structure and structure calculated according to current pocket: %1$.4f") % s_pocket_pearson
		);

		double AS_pearson = calc_correlation(
			fitting(apo_procoord, pro(POCKETAS).procoord) - apo_procoord,
			apo_pro(POCKETAS).equilibrium_coord - apo_procoord
		);
		handle_message(
			MSG_RESULT,
			boost::format("Pearson between real structure and structure calculated according to current pocket: %1$.4f") % AS_pearson
		);

		print_energy_results();
	}
	else
		handle_message(MSG_ERROR, "Lack necessary protein information.");
}

void ProAnalysis::write_matrix(MatrixXd mat, std::string writepath)
{
	std::ofstream matf(writepath, std::ios::out);
	if (matf.is_open())
	{
		matf << mat.format(CleanFmt);
		matf.close();
		handle_message(MSG_INFO, boost::format("Matrix has been written to %1%.") % writepath);
	}
	else
		handle_message(MSG_WARNING, boost::format("Can not open file %1%.") % writepath);
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
		handle_message(MSG_INFO, boost::format("Matrix (binary data) has been written to %1%.") % writepath);
	}
	else
		handle_message(MSG_WARNING, boost::format("Can not open file %1%.") % writepath);
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
		handle_message(MSG_WARNING, boost::format("Can not open file %1%.") % fpath);
}

std::string ProAnalysis::write_path(std::string fname)
{
	filesys::path workdir = workdir_path;
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

void ProAnalysis::grep_pocket_coord(VectorXd & pocket_coord, VectorXd original_coord, PocketList pocket_members)
{
	size_t i = 0;
	for (const size_t &id : pocket_members)
	{
		pocket_coord(i * 3) = original_coord(id * 3);
		pocket_coord(i * 3 + 1) = original_coord(id * 3 + 1);
		pocket_coord(i * 3 + 2) = original_coord(id * 3 + 2);
		++i;
	}
}

void ProAnalysis::modify_pocket_coord(VectorXd & coord, VectorXd replace_coord, PocketList pocket_members)
{
	size_t i = 0;
	for (const size_t &id : pocket_members)
	{
		coord(id * 3) = replace_coord(i * 3);
		coord(id * 3 + 1) = replace_coord(i * 3 + 1);
		coord(id * 3 + 2) = replace_coord(i * 3 + 2);
		++i;
	}
}

void ProAnalysis::copy_pocket_coord(VectorXd & coord, VectorXd source_coord, PocketList pocket_members)
{
	for (const size_t &id : pocket_members)
	{
		coord(id * 3) = source_coord(id * 3);
		coord(id * 3 + 1) = source_coord(id * 3 + 1);
		coord(id * 3 + 2) = source_coord(id * 3 + 2);
	}
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
		handle_message(MSG_WARNING, "Can not calculate model RMSD without pocket force generated. Call \"gen_pocket*_force\" function first.");
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
		handle_message(MSG_WARNING, "Can not calculate model correlation without pocket force generated. Call \"gen_pocket*_force\" function first.");
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

	handle_message(MSG_RESULT, buf);
}
void ProAnalysis::test_pocket(bool flag, bool info, PocketList pocket, VectorXd displacement, VectorXd pocket_force, VectorXd refcoord)
{
	if (!pocket.empty() && flag)
	{
		show_pocket(pocket);
		if (info)
		{
			handle_message(
				MSG_RESULT,
				boost::format("RMSD between real structure and structure calculated according to current pocket: %1% A.") % calc_model_rmsd(flag, pocket_force, refcoord)
			);
			handle_message(
				MSG_RESULT,
				boost::format("Pearson between real structure and structure calculated according to current pocket: %1% A.") % calc_model_correlation(flag, displacement, pocket_force)
			);
		}
		else
			handle_message(MSG_WARNING, "Lack necessary protein information.");
	}
	else
		handle_message(MSG_WARNING, "The binding pocket domain is not specificed.");
}

void ProAnalysis::show_pocket_force(bool access, PocketList pocket_members, VectorXd pocket_force)
{
	if (access)
	{
		handle_message(MSG_RESULT, "POCKET FORCE (PREDICT)");
		handle_message(MSG_EMPTY, "*---*---*---*---*---*");
		for (PocketList::const_iterator it = pocket_members.cbegin(); it != pocket_members.cend(); ++it)
		{
			Vector3d resforce = Vector3d::Zero();
			resforce << pocket_force(*it * 3), pocket_force(*it * 3 + 1), pocket_force(*it * 3 + 2);
			handle_message(MSG_EMPTY, boost::format("RES %1% FORCE %2$.4f") % (*it + 1) % calc_norm(resforce));
		}
		handle_message(MSG_EMPTY, "*---*---*---*---*---*");
	}
	else
		handle_message(MSG_WARNING, "The binding pocket domain is not specificed.");
}

void ProAnalysis::show_pro_pocket_force(PocketList pocket, VectorXd pro_force)
{
	handle_message(MSG_RESULT, "POCKET FORCE (STRUCTURE)");
	handle_message(MSG_EMPTY, "*---*---*---*---*---*");
	for (PocketList::iterator it = pocket.begin(); it != pocket.end(); ++it)
	{
		Vector3d resforce = Vector3d::Zero();
		resforce << pro_force(*it * 3), pro_force(*it * 3 + 1), pro_force(*it * 3 + 2);
		handle_message(MSG_EMPTY, boost::format("RES %1% FORCE %2$.4f") % (*it + 1) % calc_norm(resforce));
	}
	handle_message(MSG_EMPTY, "*---*---*---*---*---*");
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

	handle_message(MSG_RESULT, "PROTEIN FORCE");
	handle_message(MSG_EMPTY, "*---*---*---*---*---*");
	for (size_t i = 0; i < ProE.get_resn(); ++i)
		handle_message(MSG_EMPTY, boost::format("RES %1% FORCE %2$.4f") % (forces[i].first + 1) % forces[i].second);
	handle_message(MSG_EMPTY, "*---*---*---*---*---*");
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
			handle_message(MSG_INFO, "Given residue ID already in pocket.");
		else
		{
			pocket.push_back(id);
			handle_message(MSG_INFO, boost::format("Add residue %1% to pocket.") % (id + 1));
		}
	}
	else
		handle_message(MSG_WARNING, "Given residue ID out of range.");
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
					handle_message(MSG_INFO, boost::format("Remove residue %1% from pocket.") % (id + 1));
				}
		}
		else
			handle_message(MSG_INFO, "Given residue ID not in pocket. Can not erase it.");
	}
	else
		handle_message(MSG_WARNING, "Given residue ID out of range.");
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
			handle_message(
				MSG_WARNING,
				boost::format("Given cutoff is too short. Minimum possible cutoff is %1$.2f A.") % dist2ligand.minCoeff()
			);
	}
	else
		handle_message(MSG_WARNING, "Can not find ligand information.");
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

void ProAnalysis::calc_energy_known(bool flag, FreeEnergy & energy, PocketList pocket, VectorXd pro_force, MatrixXd distmat, VectorXd displacement)
{
	if (flag)
	{
		VectorXd force = VectorXd::Zero(pro_force.rows());
		for (PocketList::iterator it = pocket.begin(); it != pocket.end(); ++it)
		{
			force(*it * 3) = pro_force(*it * 3);
			force(*it * 3 + 1) = pro_force(*it * 3 + 1);
			force(*it * 3 + 2) = pro_force(*it * 3 + 2);
		}
		ArrayXXd distdiffmat = distmat - ProE.get_distmat();
		// Distance matrix is symmetric, so protein internal energy is 2 * 2 = 4 folds of real value
		energy.pro = (distdiffmat.pow(2) * ProE.get_kmat()).sum() / 4;
		energy.pocket = -force.transpose() * displacement;
		energy.total = energy.pro + energy.pocket;
	}
	else
		handle_message(MSG_WARNING, "Can not calculate energy without pocket force generated. Call \"gen_pocket*_force\" function first.");
}

void ProAnalysis::calc_energy_unknown(bool access, FreeEnergy & energy, VectorXd pocket_force, VectorXd equilibrium_coord)
{
	if (access)
	{
		ArrayXXd distdiffmat = gen_distmat(Dist, equilibrium_coord) - ProE.get_distmat();		
		energy.pro = (distdiffmat.pow(2) * ProE.get_kmat()).sum() / 4;
		VectorXd displacement = equilibrium_coord - ProE.get_procoord();
		energy.pocket = -pocket_force.transpose() * displacement;
		energy.total = energy.pro + energy.pocket;
	}
	else
		handle_message(MSG_WARNING, "Can not calculate energy without pocket force generated. Call \"gen_pocket*_force\" function first.");
}

void ProAnalysis::align_multiple_pockets(PocketInfo pocket1, VectorXd & pocket1_coord, PocketInfo pocket2, VectorXd & pocket2_coord)
{
	PocketList *poc1 = &(pocket1.members), *poc2 = &(pocket2.members);
	VectorXd apo_pocket1_coord = VectorXd::Zero(poc1->size() * 3), apo_pocket2_coord = VectorXd::Zero(poc2->size() * 3);
	grep_pocket_coord(apo_pocket1_coord, apo_procoord, *poc1);
	grep_pocket_coord(apo_pocket2_coord, apo_procoord, *poc2);
	pocket1_coord = fitting(apo_pocket1_coord, pocket1_coord);
	pocket2_coord = fitting(apo_pocket1_coord, pocket2_coord);

	// TODO
}

void ProAnalysis::print_energy_results()
{
	handle_message(MSG_RESULT, "FREE ENERGY RESULTS");
	handle_message(MSG_EMPTY, "*---*---*---*---*---*");
	handle_message(
		MSG_EMPTY,
		boost::format("Free energy for binding state structure S: %1$.3f Kcal/mol.") % pro(POCKETS).G.total
	);
	handle_message(
		MSG_EMPTY,
		boost::format("Pro: %1$.3f Kcal/mol.") % pro(POCKETS).G.pro
	);
	handle_message(
		MSG_EMPTY,
		boost::format("Pocket: %1$.3f Kcal/mol.") % pro(POCKETS).G.pocket
	);
	handle_message(
		MSG_EMPTY,
		boost::format("Free energy for allostery state structure A: %1$.3f Kcal/mol.") % pro(POCKETA).G.total
	);
	handle_message(
		MSG_EMPTY,
		boost::format("Pro: %1$.3f Kcal/mol.") % pro(POCKETA).G.pro
	);
	handle_message(
		MSG_EMPTY,
		boost::format("Pocket: %1$.3f Kcal/mol.") % pro(POCKETA).G.pocket
	);
	handle_message(
		MSG_EMPTY,
		boost::format("Free energy for complex state structure AS: %1$.3f Kcal/mol.") % pro(POCKETAS).G.total
	);
	handle_message(
		MSG_EMPTY,
		boost::format("Pro: %1$.3f Kcal/mol.") % pro(POCKETAS).G.pro
	);
	handle_message(
		MSG_EMPTY,
		boost::format("Pocket: %1$.3f Kcal/mol.") % pro(POCKETAS).G.pocket
	);
	handle_message(
		MSG_EMPTY,
		boost::format("Change of free energy: : %1$.3f Kcal/mol.") % ddG
	);
	handle_message(MSG_EMPTY, "*---*---*---*---*---*");
}
