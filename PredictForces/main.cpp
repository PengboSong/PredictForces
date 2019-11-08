#include <fstream>
#include <iostream>

#include <nlohmann/json.hpp>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include "defines.h"
#include "handle_io.h"
#include "Pro.h"
#include "ProAnalysis.h"

namespace filesys = boost::filesystem;
using json = nlohmann::json;

void read_settings(PredictForcesSettings &settings)
{
	json settings_json;

	std::string settings_path = "settings.json";
	if (filesys::is_regular_file(filesys::path(settings_path)))
	{
		std::ifstream file_settings(settings_path);
		if (!file_settings.is_open())
		{
			handle_message(
				MSG_ERROR,
				boost::format("Can not open settings JSON file %1.") % settings_path
			);
		}
		file_settings >> settings_json;
		file_settings.close();

		settings.workdir_path = settings_json["Workdir Path"].get<std::string>();
		settings.text_tag = settings_json["Text Tag"].get<std::string>();
		settings.binary_tag = settings_json["Binary Tag"].get<std::string>();
		settings.cache_tag = settings_json["Cache Tag"].get<std::string>();

		settings.E.pdb = settings_json["Apo State"].get<std::string>();
		settings.ES.pdb = settings_json["Binding State"].get<std::string>();
		settings.EA.pdb = settings_json["Allostery State"].get<std::string>();
		settings.EAS.pdb = settings_json["Complex State"].get<std::string>();

		for (auto& res : settings_json["Apo State Ligand Residues"])
			settings.E.ligandres.push_back(res.get<std::string>());
		for (auto& res : settings_json["Binding State Ligand Residues"])
			settings.EA.ligandres.push_back(res.get<std::string>());
		for (auto& res : settings_json["Allostery State Ligand Residues"])
			settings.ES.ligandres.push_back(res.get<std::string>());
		for (auto& res : settings_json["Complex State Ligand Residues"])
			settings.EAS.ligandres.push_back(res.get<std::string>());

		for (auto& res : settings_json["Apo State Exclude Residues"])
			settings.E.exclres.push_back(res.get<std::string>());
		for (auto& res : settings_json["Binding State Exclude Residues"])
			settings.EA.exclres.push_back(res.get<std::string>());
		for (auto& res : settings_json["Allostery State Exclude Residues"])
			settings.ES.exclres.push_back(res.get<std::string>());
		for (auto& res : settings_json["Complex State Exclude Residues"])
			settings.EAS.exclres.push_back(res.get<std::string>());

		settings.k = settings_json["k"].get<double>();
		settings.cutoff = settings_json["Cutoff"].get<double>();
	}
}

Pro init_pro(AllosteryCycleStates s, PredictForcesSettings settings)
{
	AllosteryCycleState state;
	switch (s)
	{
	case ApoState:
		state = settings.E;
		break;
	case BindingState:
		state = settings.ES;
		break;
	case AllosteryState:
		state = settings.EA;
		break;
	case ComplexState:
		state = settings.EAS;
		break;
	default:
		state = settings.E;
	}

	std::string file_pdb = state.pdb;
	if (!boost::algorithm::ends_with(file_pdb, ".pdb"))
		file_pdb += ".pdb";
	filesys::path pdb_path = filesys::path(settings.workdir_path) / file_pdb;

	ProConfigs configs = {
		settings.k,		  // Intra k
		settings.k,		  // Inter k
		settings.cutoff,  // Intra cutoff
		settings.cutoff,  // Inter cutoff
		state.ligandres,  // Ligand residues
		state.exclres     // Exclude residues
	};

	return Pro(
		pdb_path.string(),
		configs
	);
}

int main()
{
	PredictForcesSettings settings;
	read_settings(settings);
	filesys::path workdir = filesys::path(settings.workdir_path);

	if (!filesys::is_directory(workdir))
	{
		handle_message(
			MSG_ERROR,
			boost::format("Can not find directory %1%.") % workdir.string()
		);
	}

	Pro ProApo;
	Pro ProBinding;
	Pro ProAllostery;
	Pro ProComplex;

	filesys::path pdbf;

	if (settings.E.pdb.empty())
	{
		handle_message(
			MSG_ERROR,
			"Apo state structure must be loaded."
		);
	}
	else
		ProApo = init_pro(ApoState, settings);

	if (!settings.ES.pdb.empty())
		ProBinding = init_pro(BindingState, settings);

	if (!settings.EA.pdb.empty())
		ProAllostery = init_pro(AllosteryState, settings);

	if (!settings.EAS.pdb.empty())
		ProComplex = init_pro(ComplexState, settings);

	ProAnalysis Cycle;
	if (!settings.cache_tag.empty())
	{
		Cycle = ProAnalysis(
			ProApo,
			ProBinding,
			ProAllostery,
			ProComplex
		);
	}
	else
	{
		Cycle = ProAnalysis(
			ProApo,
			ProBinding,
			ProAllostery,
			ProComplex
		);
	}

	/*
	if (!settings.text_tag.empty())
	{
		Cycle.write_matrix(Cycle.get_hessian(), settings.text_tag + "_H.txt");
		Cycle.write_matrix(Cycle.get_covariance(), settings.text_tag + "_G.txt");
	}

	if (!settings.binary_tag.empty())
	{
		Cycle.write_matrix_binary(Cycle.get_hessian(), settings.binary_tag + ".hessian");
		Cycle.write_matrix_binary(Cycle.get_covariance(), settings.binary_tag + ".covariance");
	}
	*/

	Cycle.interactive();

	return 0;
}
