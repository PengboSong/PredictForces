#include "Pro.h"

Pro::Pro()
{
}

Pro::Pro(std::string fpath, bool has_ligand_flag)
{
	with_ligand_flag = has_ligand_flag;
	read(fpath);
	gen_coord();
	gen_distmat();
	gen_contact();
	pairn = contact_pairs.size();
	gen_hessian();
	gen_dist2ligand();
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

	for (size_t i = 0; i < resn; i++)
	{
		contact_map(i, i) = 1;
		for (size_t j = i + 1; j < resn; j++)
		{
			double dist_ij = 0.0;
			if (gen_distmat_flag)
				dist_ij = distmat(i, j);
			else
				dist_ij = distance(i, j);

			pair one{ i, j };
			if (pro[i].chain == pro[j].chain && distmat(i, j) < cutoff_intra)
			{
				contact_map(i, j) = contact_map(j, i) = 2;
				contact_pairs.push_back(one);
				kmat(i, j) = kmat(j, i) = k_intra;
			}
			else if (pro[i].chain != pro[j].chain && distmat(i, j) < cutoff_inter)
			{
				contact_map(i, j) = contact_map(j, i) = 3;
				contact_pairs.push_back(one);
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

void Pro::gen_hessian()
{
	hessian = Eigen::MatrixXd::Zero(3 * resn, 3 * resn);
	for (std::vector<pair>::iterator it = contact_pairs.begin(); it != contact_pairs.end(); it++)
	{
		double d = distance(*it);
		double k = k_default;
		if (get_contact(it->i, it->j) == 2)
			k = k_intra;
		else if (get_contact(it->i, it->j) == 3)
			k = k_inter;

		double hxx = -k * pow(diff_x(*it) / d, 2);
		double hyy = -k * pow(diff_y(*it) / d, 2);
		double hzz = -k * pow(diff_z(*it) / d, 2);
		double hxy = -k * diff_x(*it) * diff_y(*it) / pow(d, 2);
		double hxz = -k * diff_x(*it) * diff_z(*it) / pow(d, 2);
		double hyz = -k * diff_y(*it) * diff_z(*it) / pow(d, 2);

		hessian(3 * it->i, 3 * it->j) = hxx;
		hessian(3 * it->i, 3 * it->i) -= hxx;

		hessian(3 * it->i + 1, 3 * it->j + 1) = hyy;
		hessian(3 * it->i + 1, 3 * it->i + 1) -= hyy;

		hessian(3 * it->i + 2, 3 * it->j + 2) = hzz;
		hessian(3 * it->i + 2, 3 * it->i + 2) -= hzz;

		hessian(3 * it->i, 3 * it->j + 1) = hxy;
		hessian(3 * it->i, 3 * it->i + 1) -= hxy;
		hessian(3 * it->i + 1, 3 * it->j) = hxy;
		hessian(3 * it->i + 1, 3 * it->i) -= hxy;

		hessian(3 * it->i, 3 * it->j + 2) = hxz;
		hessian(3 * it->i, 3 * it->i + 2) -= hxz;
		hessian(3 * it->i + 2, 3 * it->j) = hxz;
		hessian(3 * it->i + 2, 3 * it->i) -= hxz;

		hessian(3 * it->i + 1, 3 * it->j + 2) = hyz;
		hessian(3 * it->i + 1, 3 * it->i + 2) -= hyz;
		hessian(3 * it->i + 2, 3 * it->j + 1) = hyz;
		hessian(3 * it->i + 2, 3 * it->i + 1) -= hyz;
	}
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

double Pro::diff_x(size_t i, size_t j)
{
	return pro[i].x - pro[j].x;
}

double Pro::diff_y(size_t i, size_t j)
{
	return pro[i].y - pro[j].y;
}

double Pro::diff_z(size_t i, size_t j)
{
	return pro[i].z - pro[j].z;
}

double Pro::distance(pair ij)
{
	return distance(ij.i, ij.j);
}

double Pro::diff_x(pair ij)
{
	return diff_x(ij.i, ij.j);
}

double Pro::diff_y(pair ij)
{
	return diff_y(ij.i, ij.j);
}

double Pro::diff_z(pair ij)
{
	return diff_z(ij.i, ij.j);
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

Eigen::MatrixXd Pro::get_hessian()
{
	return hessian;
}

Eigen::Matrix3d Pro::get_hessian(size_t i, size_t j)
{
	if (i < resn && j < resn)
		return Eigen::Matrix3d(hessian.block(3 * i, 3 * j, 3, 3));
	else
		return Eigen::Matrix3d();
}

double Pro::get_hessian_s(size_t si, size_t sj)
{
	if (si < 3 * resn && sj < 3 * resn)
		return hessian(si, sj);
	else
		return 0.0;
}

bool Pro::empty()
{
	if (resn == 0)
		return true;
	else
		return false;
}
