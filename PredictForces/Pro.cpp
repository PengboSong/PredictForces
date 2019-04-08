#include "Pro.h"

Pro::Pro()
{
}

Pro::Pro(std::string fpath, bool has_ligand_flag, std::vector<std::string> exclude, double k, double cutoff)
{
	with_ligand_flag = has_ligand_flag;
	for (std::vector<std::string>::iterator it = exclude.begin(); it != exclude.end(); ++it)
		exclres.emplace(*it);
	k_default = k_inter = k_intra = k;
	cutoff_inter = cutoff_intra = cutoff;
	std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(1);
	std::cout << "[Info] Spring constant = " << k << " J/(mol A^2)." << std::endl;
	std::cout << "[Info] Cutoff = " << cutoff << " A." << std::endl;
	read(fpath);
	std::cout << "[Info] Successfully loaded protein at path " << fpath << "." << std::endl;
	gen_coord();
	gen_distmat();
	gen_contact();
	pairn = contact_pairs.size();
	std::cout << "[Info] Coordinate matrix, distance matrix, contact map have been generated for this protein." << std::endl;
	if (with_ligand_flag)
	{
		gen_dist2ligand();
		std::cout << "[Info] Residue distance to ligand has been calculated." << std::endl;
	}
	std::cout << std::resetiosflags(std::ios::fixed);
}

Pro::~Pro()
{
}

bool Pro::has_ligand()
{
	return with_ligand_flag;
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
						std::vector<AtomInfo> grp = { read_atom(line) };
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
						std::vector<AtomInfo> grp = { read_atom(line) };
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
						std::vector<AtomInfo> grp;
						grp.push_back(ex);
						excl[resname] = grp;
					}
				}
			}
		}
		pdb.close();
		resn = proid + 1;
		proatomn = proatomid;
		ligandatomn = ligandatomid;
	}
	else
	{
		std::cout << "[Error] Unbale to open file " << fpath << '.' << std::endl;
	}
}

void Pro::gen_contact()
{
	contact_map = Eigen::MatrixXi::Zero(resn, resn);
	kmat = Eigen::ArrayXXd::Zero(resn, resn);

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
			{
				contact_map(i, j) = contact_map(j, i) = 0;
			}
		}
	}

	gen_contact_flag = true;
}

void Pro::gen_coord()
{
	procoord = Eigen::VectorXd::Zero(3 * resn);
	for (size_t i = 0; i < resn; i++)
	{
		procoord(3 * i) = pro[i].x;
		procoord(3 * i + 1) = pro[i].y;
		procoord(3 * i + 2) = pro[i].z;
	}

	ligandcoord = Eigen::VectorXd::Zero(3 * ligandatomn);
	size_t j = 0;
	for (std::map<std::string, std::vector<AtomInfo>>::iterator it = ligand.begin(); it != ligand.end(); ++it)
	{
		for (std::vector<AtomInfo>::iterator iit = it->second.begin(); iit != it->second.end(); ++iit)
		{
			ligandcoord(3 * j) = iit->x;
			ligandcoord(3 * j + 1) = iit->y;
			ligandcoord(3 * j + 2) = iit->z;
			++j;
		}
	}

	for (std::map<size_t, std::vector<AtomInfo>>::iterator it = proatoms.begin(); it != proatoms.end(); ++it)
	{
		Eigen::VectorXd rescoord = Eigen::VectorXd::Zero(3 * it->second.size());
		size_t k = 0;
		for (std::vector<AtomInfo>::iterator iit = it->second.begin(); iit != it->second.end(); ++iit)
		{
			rescoord(3 * k) = iit->x;
			rescoord(3 * k + 1) = iit->y;
			rescoord(3 * k + 2) = iit->z;
			++k;
		}
		rescoords[it->first] = rescoord;
	}
}

Eigen::MatrixXd Pro::gen_hessian()
{
	Eigen::MatrixXd hessian = Eigen::MatrixXd::Zero(3 * resn, 3 * resn);

	if (gen_contact_flag)
	{
		for (std::vector<std::pair<size_t, size_t>>::iterator it = contact_pairs.begin(); it != contact_pairs.end(); ++it)
		{
			size_t pi = it->first;
			size_t pj = it->second;
			double diffx = pro[pi].x - pro[pj].x;
			double diffy = pro[pi].y - pro[pj].y;
			double diffz = pro[pi].z - pro[pj].z;
			double d = pow(diffx, 2) + pow(diffy, 2) + pow(diffz, 2);
			double k = k_default;
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

Eigen::MatrixXd Pro::gen_covariance(Eigen::MatrixXd hessian)
{
	Eigen::MatrixXd covariance = Eigen::MatrixXd::Zero(3 * resn, 3 * resn);

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(hessian);
	Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();
	Eigen::VectorXd zero2inf_eigenvalues = eigenvalues;
	Eigen::MatrixXd eigenvectors = eigensolver.eigenvectors();
	std::vector<size_t> zeromodes, nonzeromodes;
	size_t zeromoden = calc_zero_modes(eigenvalues, zero2inf_eigenvalues);

	if (zeromoden == 6)
	{
		Eigen::MatrixXd U = Eigen::ArrayXXd(eigenvectors).colwise() / Eigen::ArrayXd(zero2inf_eigenvalues);
		covariance = (kB * Navo * Temp / k_default) * (eigenvectors.transpose() * U);
		/*
		for (size_t i = 0; i < 3 * resn; ++i)
			for (size_t j = 0; j < 3 * resn; ++j)
				for (std::vector<size_t>::iterator k = nonzeromodes.begin(); k != nonzeromodes.end(); ++k)
					covariance(i, j) += eigenvectors(*k, i) * eigenvectors(*k, j) / eigenvalues(*k);
		*/
	}
	else
		std::cout << "[Error] Hessian matrix has " << zeromoden << " zero modes. Please check it before constructing covariance matrix." << std::endl;

	return covariance;
}

void Pro::gen_distmat()
{
	distmat = Eigen::MatrixXd::Zero(resn, resn);
	for (size_t i = 0; i < resn; i++)
		for (size_t j = i + 1; j < resn; j++)
			distmat(j, i) = distmat(i, j) = distance(i, j);

	gen_distmat_flag = true;
}

void Pro::gen_dist2ligand()
{
	dist2ligand = Eigen::VectorXd::Zero(resn);

	if (with_ligand_flag)
	{
		Eigen::Map<Eigen::Matrix3Xd> ligandxyz(ligandcoord.data(), 3, ligandatomn);
		for (size_t i = 0; i < resn; ++i)
		{
			size_t resatomn = get_resatomn(i);
			Eigen::VectorXd rescoord = get_rescoord(i);
			Eigen::Map<Eigen::Matrix3Xd> resxyz(rescoord.data(), 3, resatomn);
			Eigen::VectorXd dist = Eigen::VectorXd::Zero(resatomn);
			Eigen::Array3Xd diffxyz(3, ligandxyz.cols());
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

size_t Pro::calc_zero_modes(Eigen::VectorXd eigenvalues, std::vector<size_t> *zeromodes, std::vector<size_t> *nonzeromodes)
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

size_t Pro::calc_zero_modes(Eigen::VectorXd eigenvalues, Eigen::VectorXd &zero2inf_eigenvalues)
{
	size_t count = 0;

	zero2inf_eigenvalues = eigenvalues;

	for (size_t i = 0; i < size_t(eigenvalues.size()); i++)
		if (eigenvalues(i) == 0.0 || abs(eigenvalues(i)) < 1e-10)
		{
			zero2inf_eigenvalues(i) = std::numeric_limits<double>::infinity();
			++count;
		}
	return count;
}

bool Pro::has_res(size_t id)
{
	if (id < resn)
		return true;
	else
		return false;
}

std::string Pro::get_resname(size_t id)
{
	if (id < resn)
		return pro[id].resname;
	else
		return std::string();
}

std::string Pro::get_chain(size_t id)
{
	if (id < resn)
		return pro[id].chain;
	else
		return std::string();
}

size_t Pro::get_resid(size_t id)
{
	if (id < resn)
		return pro[id].resid;
	else
		return 0;
}

double Pro::get_x(size_t id)
{
	if (id < resn)
		return pro[id].x;
	else
		return 0.0;
}

double Pro::get_y(size_t id)
{
	if (id < resn)
		return pro[id].y;
	else
		return 0.0;
}

double Pro::get_z(size_t id)
{
	if (id < resn)
		return pro[id].z;
	else
		return 0.0;
}

double Pro::get_bfactor(size_t id)
{
	if (id < resn)
		return pro[id].bfactor;
	else
		return 0.0;
}

size_t Pro::get_resn()
{
	return resn;
}

int Pro::get_contact(size_t i, size_t j)
{
	if (i < resn && j < resn)
		return contact_map(i, j);
	else
		return 0;
}

void Pro::show_contact_pairs()
{
	for (std::vector<std::pair<size_t, size_t>>::iterator it = contact_pairs.begin(); it != contact_pairs.end(); ++it)
	{
		std::cout << '(' << it->first << ',' << it->second << ')' << std::endl;
	}
}

Eigen::MatrixXi Pro::get_contact_map()
{
	return contact_map;
}

Eigen::VectorXd Pro::get_procoord()
{
	return procoord;
}

Eigen::VectorXd Pro::get_ligandcoord()
{
	return ligandcoord;
}

Eigen::VectorXd Pro::get_dist2ligand()
{
	return dist2ligand;
}

Eigen::VectorXd Pro::get_rescoord(size_t id)
{
	if (id < resn)
		return rescoords[id];
	else
		return Eigen::VectorXd();
}

size_t Pro::get_resatomn(size_t id)
{
	if (id < resn)
		return rescoords[id].size() / 3;
	else
		return 0;
}

Eigen::MatrixXd Pro::get_distmat()
{
	return distmat;
}

Eigen::ArrayXXd Pro::get_kmat()
{
	return kmat;
}

bool Pro::empty()
{
	if (resn == 0)
		return true;
	else
		return false;
}
