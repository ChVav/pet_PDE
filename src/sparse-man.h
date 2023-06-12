#pragma once

#include "common.h"
#include <vector>
#include <string>

class sparse_mat { // manual implementation for sparse matrix

public:
	std::vector<double> vector_mult(std::vector<double>, double);

	std::vector<int> _Get_pos() {
		return _pos;
	}
	std::vector<double> _Get_val() {
		return _val;
	}

	sparse_mat();
	sparse_mat(int dim, std::vector<int> pos, std::vector<double> val);

private:
	int _dim;
	int _count;
	std::vector<int> _pos;
	std::vector<double> _val;
};

sparse_mat assemble_matrix(
	const std::vector<point>& points,
	const std::vector<triangle>& triangles,
	const std::vector<line>& lines,
	const std::vector<double>& areas);

std::vector<double> one_timestep(
	double dt,
	const sparse_mat& B,
	const std::vector<double>& u_n);