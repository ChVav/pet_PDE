
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

Eigen::VectorXd one_timestep(double dt, Eigen::SparseMatrix<double>& B, Eigen::VectorXd u_n) {
	auto u_n1 = u_n + dt * (B * u_n);
	return u_n1;
}

int main() {
	double dt = 0.01;

	Eigen::SparseMatrix<double> B(4, 4);
	B.insert(0, 0) = 1.0;
	B.insert(0, 3) = 4.0;
	B.insert(1, 1) = 3.0;
	B.insert(3, 0) = 5.0;

	std::cout << "sparse matrix B = " << B << std::endl;

	std::vector<double> ini_state = { 0.1,0.1,0.1,0.1 };

	std::cout << "ini_state first element = " << ini_state[0] << std::endl;

	Eigen::VectorXd u = Eigen::Map<Eigen::VectorXd>(ini_state.data(), ini_state.size());

	std::cout << "u = " << u << std::endl;

	Eigen::VectorXd u1 = u + dt * B * u;

	std::cout << "u1 = " << u1 << std::endl;

	Eigen::VectorXd u2 = one_timestep(dt, B, u1);

	std::cout << "u2 = " << u2 << std::endl;

	return 0;
}