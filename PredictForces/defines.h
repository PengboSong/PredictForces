#ifndef PREDICT_FORCES_DEFINES_H
#define PREDICT_FORCES_DEFINES_H

#include <string>
#include <vector>
#include <list>
#include <map>
#include <utility>

#include <boost/format.hpp>
#include <boost/filesystem.hpp>
namespace filesys = boost::filesystem;

#include <Eigen/Dense>
using namespace Eigen;

#include "HandleMessage.h"

/* <--- Constants ---> */

constexpr double Navo = 6.02214076e23;
constexpr double kB = 1.380649e-23;
constexpr double Temp = 298.15;
constexpr double h = 6.6260693e-34;
constexpr double PI = 3.1415926535897932;

/* <--- Enumeration ---> */

enum AllosteryCycleStates : uint8_t {
	ApoState,		 // E
	BindingState,	 // ES
	AllosteryState,	 // EA
	ComplexState	 // EAS
};

enum Pockets : uint8_t {
	POCKETS,	  // Binding Pocket
	POCKETA,      // Allostery Pocket
	POCKETAS	  // Complex Pocket
};

enum LFmethods : uint8_t {
	BatchGradientDescent,
	NormalEquation
};

enum DistMatType : uint8_t {
	Dist,
	XYZdiff,
	XYZdiffLtCutoff
};

enum EMType : uint8_t {
	FixPocket,
	FixPocketWithForce,
	OnlyForce
};

/* <--- Alias ---> */

typedef std::list<size_t> PocketList;

/* <--- Struct ---> */

typedef struct {
	double learning_rate = 0.0;
	double convergence = 0.0;
	double upward_factor = 0.0;
	double downward_factor = 0.0;
	size_t niteration = 0;
} BGDpara;

typedef struct
{
	std::string resname;
	std::string chain;
	size_t resid;
	double x;
	double y;
	double z;
	double bfactor;
} ResInfo;

typedef struct
{
	std::string atomname;
	std::string resname;
	std::string chain;
	size_t atomid;
	size_t resid;
	double x;
	double y;
	double z;
	double bfactor;
} AtomInfo;

typedef struct
{
	std::string pdb;
	std::vector<std::string> ligandres;
	std::vector<std::string> exclres;
} AllosteryCycleState;

typedef struct
{
	std::string workdir_path;
	std::string text_tag;
	std::string binary_tag;
	std::string cache_tag;
	AllosteryCycleState E;
	AllosteryCycleState ES;
	AllosteryCycleState EA;
	AllosteryCycleState EAS;
	double k;
	double cutoff;
} PredictForcesSettings;

typedef struct
{
	double k_intra;
	double k_inter;
	double cutoff_intra;
	double cutoff_inter;
	std::vector<std::string> ligandres;
	std::vector<std::string> exclres;
} ProConfigs;

typedef struct {
	bool access_force = false;
	double rmsd = 0.0;
	PocketList members;
	VectorXd force;
	VectorXd coord;
	MatrixXd displacement;
} PocketInfo;

typedef struct {
	double pro = 0.0;
	double pocket = 0.0;
	double total = 0.0;
} FreeEnergy;

typedef struct {
	bool empty = true;
	bool withligand = false;
	bool preprocess = false;
	double rmsd = 0.0;
	double meanforce = 0.0;
	double pearson = 0.0;
	VectorXd coord;
	VectorXd fitcoord;
	VectorXd dist2ligand;
	VectorXd equilibrium_coord;
	VectorXd displacement;
	VectorXd force;
	MatrixXd distdiff;
	FreeEnergy G;
} ProInfo;

/* <--- Alias ---> */

typedef std::vector<AtomInfo> AtomInfoList;

typedef std::pair<size_t, size_t> ContactPair;

typedef std::map<Pockets, PocketInfo> PocketContainer;
typedef std::map<Pockets, ProInfo> ProInfoContainer;

#endif
