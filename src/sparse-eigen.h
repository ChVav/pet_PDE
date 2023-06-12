#pragma once

#include <Eigen/Sparse>
#include "common.h"
#include <vector>
#include <string>

Eigen::SparseMatrix<double> assemble_matrix(
	const std::vector<point>& points,
	const std::vector<triangle>& triangles,
	const std::vector<line>& lines,
	const std::vector<double>& areas);

Eigen::VectorXd one_timestep(
	double dt,
	const Eigen::SparseMatrix<double>& B,
	const Eigen::VectorXd& u_n);

//overload output function so that also values from Eigen::VectorXd can be written to output
void output(
	double time_step,
	const std::vector<point>& points,
	const std::vector<line>& lines,
	const std::vector<triangle>& triangles,
	Eigen::VectorXd& u);