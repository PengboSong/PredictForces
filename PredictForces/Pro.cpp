#include "Pro.h"

Pro::Pro()
{
}

Pro::Pro(std::string fpath, ProConfigs configs)
{
	if (configs.ligandres.size() != 0)
	{
		with_ligand_flag = true;
		for (std::vector<std::string>::iterator it = configs.ligandres.begin(); it != configs.ligandres.end(); ++it)
			ligandres.emplace(*it);
	}
	if (configs.exclres.size() != 0)
	{
		for (std::vector<std::string>::iterator it = configs.exclres.begin(); it != configs.exclres.end(); ++it)
			exclres.emplace(*it);
	}

	k_intra = configs.k_intra;
	k_inter = configs.k_inter;
	if (k_intra == k_inter)
	{
		handle_message(
			MSG_INFO,
			boost::format("Spring constant = %1$.1f Kcal/(mol A^2).") % k_intra
		);
	}
	else
	{
		handle_message(
			MSG_INFO,
			boost::format("Intra Spring constant = %1$.1f Kcal/(mol A^2).") % k_intra
		);
		handle_message(
			MSG_INFO,
			boost::format("Inter Spring constant = %1$.1f Kcal/(mol A^2).") % k_inter
		);
	}

	cutoff_intra = configs.cutoff_intra;
	cutoff_inter = configs.cutoff_inter;
	if (cutoff_intra == cutoff_inter)
	{
		handle_message(
			MSG_INFO,
			boost::format("Cutoff = %1$.1f A.") % cutoff_intra
		);
	}
	else
	{
		handle_message(
			MSG_INFO,
			boost::format("Intra Cutoff = %1$.1f A.") % cutoff_intra
		);
		handle_message(
			MSG_INFO,
			boost::format("Inter Cutoff = %1$.1f A.") % cutoff_inter
		);
	}

	read(fpath);
	handle_message(
		MSG_INFO,
		boost::format("Successfully loaded protein at path %1%") % fpath
	);

	gen_coord();
	gen_distmat();
	gen_contact();
	pairn = contact_pairs.size();
	handle_message(
		MSG_INFO,
		"Coordinate matrix, distance matrix, contact map have been generated for this protein."
	);

	if (with_ligand_flag)
	{
		gen_dist2ligand();
		handle_message(
			MSG_INFO,
			"Residue distance to ligand has been calculated."
		);
	}
}

Pro::~Pro()
{
}

void Pro::read(std::string fpath)
{
	std::string line;
	std::ifstream pdb(fpath);
	size_t proid = 0, proatomid = 0, ligandatomid = 0;
	std::string prev_chain = "";
	size_t prev_resid = 0;
	if (pdb.is_open())
	{
		while (std::getline(pdb, line))
		{
			std::string record = read_record(line);
			if (record == "ATOM" || record == "HETATM")
			{
				std::string resname = read_resname(line);

				if (prores.find(resname) != prores.end())
				{
					if (prev_resid == 0 && prev_chain == "")
					{
						prev_resid = read_resid(line);
						prev_chain = read_chain(line);
					}

					if (prev_resid != read_resid(line) || prev_chain != read_chain(line))
						++proid;

					if (read_atomname(line) == "CA")
					{
						ResInfo res = read_res(line);
						pro[proid] = res;
					}

					if (proatoms.find(proid) != proatoms.end())
						proatoms[proid].push_back(read_atom(line));
					else
					{
						AtomInfoList grp = { read_atom(line) };
						proatoms[proid] = grp;
					}

					prev_resid = read_resid(line);
					prev_chain = read_chain(line);
					++proatomid;
				}
				else if (exclres.find(resname) == exclres.end())
				{
					if (ligand.find(resname) != ligand.end())
						ligand[resname].push_back(read_atom(line));
					else
					{
						AtomInfoList grp = { read_atom(line) };
						ligand[resname] = grp;
					}

					++ligandatomid;
				}
				else
				{
					AtomInfo ex = read_atom(line);

					if (excl.find(resname) != excl.end())
						excl[resname].push_back(ex);
					else
					{
						AtomInfoList grp;
						grp.push_back(ex);
						excl[resname] = grp;
					}
				}
			}
		}
		pdb.close();
		line.clear();
		resn = proid + 1;
		proatomn = proatomid;
		ligandatomn = ligandatomid;
	}
	else
		handle_message(
			MSG_ERROR,
			boost::format("Unbale to open file %1%.") % fpath
		);
}

void Pro::gen_contact()
{
	contact_map = MatrixXi::Zero(resn, resn);
	kmat = ArrayXXd::Zero(resn, resn);

	for (size_t i = 0; i < resn; ++i)
	{
		contact_map(i, i) = 1;
		for (size_t j = i + 1; j < resn; ++j)
		{
			double dist_ij = gen_distmat_flag ? distmat(i, j) : distance(i, j);

			if (pro[i].chain == pro[j].chain && dist_ij < cutoff_intra)
			{
				contact_map(i, j) = contact_map(j, i) = 2;
				contact_pairs.push_back(std::make_pair(i, j));
				contact_pairs.push_back(std::make_pair(j, i));
				kmat(i, j) = kmat(j, i) = k_intra;
			}
			else if (pro[i].chain != pro[j].chain && dist_ij < cutoff_inter)
			{
				contact_map(i, j) = contact_map(j, i) = 3;
				contact_pairs.push_back(std::make_pair(i, j));
				contact_pairs.push_back(std::make_pair(j, i));
				kmat(i, j) = kmat(j, i) = k_inter;
			}
			else
				contact_map(i, j) = contact_map(j, i) = 0;
		}
	}

	gen_contact_flag = true;
}

void Pro::gen_coord()
{
	procoord = VectorXd::Zero(3 * resn);
	for (size_t i = 0; i < resn; i++)
	{
		procoord(3 * i) = pro[i].x;
		procoord(3 * i + 1) = pro[i].y;
		procoord(3 * i + 2) = pro[i].z;
	}

	ligandcoord = VectorXd::Zero(3 * ligandatomn);
	size_t j = 0;
	for (std::map<std::string, AtomInfoList>::iterator it = ligand.begin(); it != ligand.end(); ++it)
	{
		for (AtomInfoList::iterator iit = it->second.begin(); iit != it->second.end(); ++iit)
		{
			ligandcoord(3 * j) = iit->x;
			ligandcoord(3 * j + 1) = iit->y;
			ligandcoord(3 * j + 2) = iit->z;
			++j;
		}
	}

	for (std::map<size_t, AtomInfoList>::iterator it = proatoms.begin(); it != proatoms.end(); ++it)
	{
		VectorXd rescoord = VectorXd::Zero(3 * it->second.size());
		size_t k = 0;
		for (AtomInfoList::iterator iit = it->second.begin(); iit != it->second.end(); ++iit)
		{
			rescoord(3 * k) = iit->x;
			rescoord(3 * k + 1) = iit->y;
			rescoord(3 * k + 2) = iit->z;
			++k;
		}
		rescoords[it->first] = rescoord;
	}
}

MatrixXd Pro::gen_hessian()
{
	MatrixXd hessian = MatrixXd::Zero(3 * resn, 3 * resn);

	if (gen_contact_flag)
	{
		for (std::vector<ContactPair>::iterator it = contact_pairs.begin(); it != contact_pairs.end(); ++it)
		{
			size_t pi = it->first;
			size_t pj = it->second;
			double diffx = pro[pi].x - pro[pj].x;
			double diffy = pro[pi].y - pro[pj].y;
			double diffz = pro[pi].z - pro[pj].z;
			double d = pow(diffx, 2) + pow(diffy, 2) + pow(diffz, 2);

			double k = k_intra;
			if (get_contact(pi, pj) == 2)
				k = k_intra;
			else if (get_contact(pi, pj) == 3)
				k = k_inter;
			k /= d;

			double hxx = -k * pow(diffx, 2);
			double hyy = -k * pow(diffy, 2);
			double hzz = -k * pow(diffz, 2);
			double hxy = -k * diffx * diffy;
			double hxz = -k * diffx * diffz;
			double hyz = -k * diffy * diffz;

			hessian(3 * pi, 3 * pj) = hxx;
			hessian(3 * pi, 3 * pi) -= hxx;

			hessian(3 * pi + 1, 3 * pj + 1) = hyy;
			hessian(3 * pi + 1, 3 * pi + 1) -= hyy;

			hessian(3 * pi + 2, 3 * pj + 2) = hzz;
			hessian(3 * pi + 2, 3 * pi + 2) -= hzz;

			hessian(3 * pi, 3 * pj + 1) = hxy;
			hessian(3 * pi + 1, 3 * pj) = hxy;
			hessian(3 * pi, 3 * pi + 1) -= hxy;
			hessian(3 * pi + 1, 3 * pi) -= hxy;

			hessian(3 * pi, 3 * pj + 2) = hxz;
			hessian(3 * pi + 2, 3 * pj) = hxz;
			hessian(3 * pi, 3 * pi + 2) -= hxz;
			hessian(3 * pi + 2, 3 * pi) -= hxz;

			hessian(3 * pi + 1, 3 * pj + 2) = hyz;
			hessian(3 * pi + 2, 3 * pj + 1) = hyz;
			hessian(3 * pi + 1, 3 * pi + 2) -= hyz;
			hessian(3 * pi + 2, 3 * pi + 1) -= hyz;
		}
	}
	return hessian;
}

MatrixXd Pro::gen_covariance(MatrixXd hessian)
{
	MatrixXd covariance = MatrixXd::Zero(3 * resn, 3 * resn);

	SelfAdjointEigenSolver<MatrixXd> eigensolver(hessian);
	VectorXd eigenvalues = eigensolver.eigenvalues();
	VectorXd zero2inf_eigenvalues = eigenvalues;
	MatrixXd eigenvectors = eigensolver.eigenvectors();
	std::vector<size_t> zeromodes, nonzeromodes;
	size_t zeromoden = calc_zero_modes(eigenvalues, zero2inf_eigenvalues);

	if (zeromoden == 6)
	{
		MatrixXd U = ArrayXXd(eigenvectors.transpose()).colwise() / ArrayXd(zero2inf_eigenvalues);
		covariance = (kB * Navo * Temp / k_intra) * (eigenvectors * U);
				
		/*
		for (size_t i = 0; i < 3 * resn; ++i)
			for (size_t j = 0; j < 3 * resn; ++j)
				for (vector<size_t>::iterator k = nonzeromodes.begin(); k != nonzeromodes.end(); ++k)
					covariance(i, j) += eigenvectors(*k, i) * eigenvectors(*k, j) / eigenvalues(*k);
		
		for (size_t i = 0; i < 3 * resn; ++i)
			for (size_t j = 0; j <= i; ++j)
				for (size_t k = 0; k < 3 * resn; ++k)
					if (eigenvalues(k) > 1e-5)
					{
						covariance(i, j) += (kB * Navo * Temp / k_default) * eigenvectors(k, i) * eigenvectors(k, j) / eigenvalues(k);
						covariance(j, i) = covariance(i, j);
					}
		*/
	}
	else
		handle_message(
			MSG_ERROR,
			boost::format("Hessian matrix has %1% zero modes. Please check it before constructing covariance matrix.") % zeromoden
		);

	return covariance;
}

void Pro::gen_entropy(MatrixXd hessian, double &entropy)
{
	SelfAdjointEigenSolver<MatrixXd> eigensolver(hessian);
	VectorXd eigenvalues = eigensolver.eigenvalues();
	VectorXd zero2inf_eigenvalues = eigenvalues;
	std::vector<size_t> zeromodes, nonzeromodes;
	size_t zeromoden = calc_zero_modes(eigenvalues, zero2inf_eigenvalues);
	for (size_t i = 0; i < eigenvalues.size(); ++i)
		entropy -= kB * Temp * log(kB * Temp / h / sqrt(zero2inf_eigenvalues(i)));
}

void Pro::gen_distmat()
{
	distmat = MatrixXd::Zero(resn, resn);
	for (size_t i = 0; i < resn; i++)
		for (size_t j = i + 1; j < resn; j++)
			distmat(j, i) = distmat(i, j) = distance(i, j);

	gen_distmat_flag = true;
}

void Pro::gen_dist2ligand()
{
	dist2ligand = VectorXd::Zero(resn);

	if (with_ligand_flag)
	{
		Map<Matrix3Xd> ligandxyz(ligandcoord.data(), 3, ligandatomn);
		for (size_t i = 0; i < resn; ++i)
		{
			size_t resatomn = get_resatomn(i);
			VectorXd rescoord = get_rescoord(i);
			Map<Matrix3Xd> resxyz(rescoord.data(), 3, resatomn);
			VectorXd dist = VectorXd::Zero(resatomn);
			Array3Xd diffxyz(3, ligandxyz.cols());
			for (size_t j = 0; j < resatomn; ++j)
			{
				diffxyz = ligandxyz.colwise() - resxyz.col(j);
				dist(j) = sqrt(diffxyz.pow(2).colwise().sum().minCoeff());
			}
			dist2ligand(i) = dist.minCoeff();
		}
	}
}

double Pro::distance(size_t i, size_t j)
{
	return sqrt(pow(pro[i].x - pro[j].x, 2) + pow(pro[i].y - pro[j].y, 2) + pow(pro[i].z - pro[j].z, 2));
}

size_t Pro::calc_zero_modes(VectorXd eigenvalues, std::vector<size_t> *zeromodes, std::vector<size_t> *nonzeromodes)
{
	size_t count = 0;
	for (size_t i = 0; i < size_t(eigenvalues.size()); i++)
	{
		if (eigenvalues(i) == 0.0 || abs(eigenvalues(i)) < 1e-10 )
		{
			++count;
			zeromodes->push_back(i);
		}
		else
		{
			nonzeromodes->push_back(i);
		}
	}
	return count;
}

size_t Pro::calc_zero_modes(VectorXd eigenvalues, VectorXd &zero2inf_eigenvalues)
{
	size_t count = 0;

	zero2inf_eigenvalues = eigenvalues;

	for (size_t i = 0; i < size_t(eigenvalues.size()); i++)
		if (eigenvalues(i) == 0.0 || abs(eigenvalues(i)) < 1e-8)
		{
			zero2inf_eigenvalues(i) = std::numeric_limits<double>::infinity();
			++count;
		}
	return count;
}

void Pro::show_contact_pairs()
{
	std::string pairsbuf;
	boost::format pairsformat("(%1%, %2%)");
	for (std::vector<ContactPair>::iterator it = contact_pairs.begin(); it != contact_pairs.end(); ++it)
		pairsbuf += (pairsformat % it->first % it->second).str();
	handle_message(
		MSG_RESULT,
		"Contact pairs:" + pairsbuf
	);
}
