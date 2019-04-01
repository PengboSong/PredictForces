#include <filesystem>

#include "Pro.h"
#include "ProAnalysis.h"

std::list<double> cutoff_pocket = { 2.5, 3.0, 3.5, 4.0, 4.5, 5.0 };

void write_binary(const char* filename, const Eigen::MatrixXd &matrix) {
	std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
	typename Eigen::MatrixXd::Index rows = matrix.rows(), cols = matrix.cols();
	out.write((char*)(&rows), sizeof(typename Eigen::MatrixXd::Index));
	out.write((char*)(&cols), sizeof(typename Eigen::MatrixXd::Index));
	out.write((char*)matrix.data(), rows*cols * sizeof(typename Eigen::MatrixXd::Scalar));
	out.close();
}

int main(int argc, char** argv)
{
	std::string tmp;
	std::string dataset = "C:\\Users\\china\\source\\repos\\PredictForces\\dataset\\";
	std::string metjapo = dataset + "MetJ\\1cmb.pdb";
	std::string metjc = dataset + "MetJ\\1cma.pdb";
	std::string metja = dataset + "MetJ\\1cmc.pdb";

	std::string p38apo = dataset + "p38\\1r39.pdb";
	std::string p38a1 = dataset + "p38\\1kv1.pdb";
	std::string p38a2 = dataset + "p38\\1kv2.pdb";

	Pro Apo(metjapo, false);
	Pro Alle(metja, true);
	Pro Complex(metjc, true);

	ProAnalysis Metj(Apo, Pro::Pro(), Alle, Complex);

	Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
	std::ofstream hessianf("hessian.txt");
	hessianf << Metj.get_hessian_matrix().format(CleanFmt);
	std::ofstream convariancef("convariance.txt");
	convariancef << Metj.get_convariance_matrix().format(CleanFmt);
	hessianf.close();
	convariancef.close();

	Metj.set_learning_step(1e-6);
	Metj.set_convergence(1e-6);
	Metj.set_iteration_times(1000000);
	Metj.gen_pocketA(3.5);
	Metj.test_pocketA();

	Metj.set_learning_step(1e-6);
	Metj.set_convergence(1e-6);
	Metj.set_iteration_times(1000000);
	Metj.gen_pocketAS(3.5);
	Metj.test_pocketAS();
	
	std::cin >> tmp;
	return 0;
}