#pragma once

#include "common.h"
#include <vector>
#include <string>

class sparse_mat { // manual csr-implementation for sparse matrix

public:
	std::vector<double> vector_mult(std::vector<double>);

	std::vector<int> _Get_col_ind() {
        return _col_ind;
    }
	std::vector<int> _Get_row_ind() {
        return _row_ind;
    }
	std::vector<int> _Get_row_ptr() {
        return _row_ptr;
    }
    std::vector<double> _Get_val() {
        return _val;
	}
	sparse_mat();
    sparse_mat(int, std::vector<int>, std::vector<int>, std::vector<double>);
private:
	int _dim;
    int _count;
    std::vector<int> _col_ind;
    std::vector<int> _row_ind;
	std::vector<int> _row_ptr;
    std::vector<double> _val;
};

std::vector<int> convert_coo_to_csr(int, std::vector<int>);

sparse_mat assemble_matrix(
	const std::vector<point>& points,
	const std::vector<triangle>& triangles,
	const std::vector<line>& lines,
	const std::vector<double>& areas);

std::vector<double> one_timestep(
	double dt,
	const sparse_mat& B,
	const std::vector<double>& u_n);