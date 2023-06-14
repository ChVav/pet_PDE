
/* Benchmarking sparse-csr implementation*/

#include <benchmark/benchmark.h>
#include "common.h"
#include "sparse-csr.h"

#include <cassert>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>

using namespace std;

//Functions unique to sparse-csr implementation
vector<int> convert_coo_to_csr(int dim, vector<int> row_ind) {
	vector<int> row_ptr(dim+1, 0);
	for (int i = 0; i < row_ind.size(); i++) {
    	row_ptr[row_ind[i] + 1]++;
	}
	for (int i = 0; i < dim; i++) {
    	row_ptr[i + 1] += row_ptr[i];
	}
	return row_ptr;
}

sparse_mat::sparse_mat() {}
sparse_mat::sparse_mat(int dim, vector<int> col_ind, vector<int> row_ind, vector<double> val) {
        _count=col_ind.size();
		_dim=dim;
		_col_ind=col_ind;
		_row_ind=row_ind;
		_val=val;
		_row_ptr=convert_coo_to_csr(dim, row_ind);
    }

vector<double> sparse_mat::vector_mult(vector<double> u) {
    vector<double> res(u.size());
    for (int i=0; i<u.size(); i++) {
		res[i]=0;
		for (int j= _row_ptr[i]; j<_row_ptr[i+1]; j++) {
			res[i]+= _val[j]*u[_col_ind[j]];
		}
	}
    return res;
}

// calculate matrix B
sparse_mat assemble_matrix(
	const vector<point>& points,
	const vector<triangle>& triangles,
	const vector<line>& lines,
	const vector<double>& areas) {

	int n_v = points.size();
	int n_T = triangles.size();


	// use pseudo code from finite_elements.pdf:
	vector<int> col_ind;
	vector<int> row_ind;
	vector<double> values;
	for (int i = 0; i < n_v; i++) {
		double H = 0;
		for (int k = 0; k < n_T; k++) {
			if (triangles[k].has_vertex(i)) {
				H += areas[k];
			}
		}

		if (!on_boundary(i, lines)) {
			for (int j = 0; j < n_v; j++) {
				double Bij = 0;
				for (int k = 0; k < n_T; k++) {
					if (triangles[k].has_vertex(i) && triangles[k].has_vertex(j)) {
						double bik, bjk;
						double cik, cjk;
						compute_bc(i, triangles[k], points, bik, cik);
						compute_bc(j, triangles[k], points, bjk, cjk);
						Bij -= areas[k] * (bik * bjk + cik * cjk);
					}
				}
				if (Bij != 0) {
					Bij /= H;
					row_ind.push_back(i);
					col_ind.push_back(j);
					values.push_back(Bij);
				}
			}
		}
	}
	// construct B
	sparse_mat B;
	B = sparse_mat(points.size(), col_ind, row_ind, values);
	return B;
}

// forward Euler, calculate one timestep
vector<double> one_timestep(double dt, sparse_mat& B, const vector<double>& u_n) {
	//returns u_n+1=u_n+dt*B*u_n
	vector<double> u_n1(u_n.size());
	vector<double> delta_u(u_n.size());
	delta_u=B.vector_mult(u_n);
	for (int i=0; i<u_n.size(); i++) {
		u_n1[i]=dt*delta_u[i]+u_n[i];
	}
	return u_n1;
}

// Combined functions for benchmarking sparse-man implementation (assemble_matrix + forward Euler + initial heat state): 
void pdeSparseCsr(const vector<point>& points,
	const vector<triangle>& triangles,
	const vector<line>& lines,
	const vector<double>& areas,
	double dt) {

	sparse_mat B;
	B = assemble_matrix(points, triangles, lines, areas);

	// actual time evolution:
	//initial state
	vector<double> ini_state(points.size());
	for (int i = 0; i < points.size(); i++) {
		double x = points[i].Getx(), y = points[i].Gety();
		ini_state[i] = heat(x, y);
	}

	vector<double> u = ini_state;
	int count = 0;
	int file_count = 1;
	double end_time = 0.1;
	int frame_count = int(end_time / dt);
	for (int i = 0; i < frame_count; i++) {
		//double dt=0.01;
		u = one_timestep(dt, B, u); //calculates one timestep
		count += 1;
	}
}

// testing assembly of sparse-man matrix only
// fixed mesh
static void assembleMatrixSparseCsrBench1(benchmark::State& s) {
	// Input data mesh
	vector<point> points;
	vector<triangle> triangles;
	vector<line> lines;
	read_mesh("assets/square.msh", points, triangles, lines);

	//Calculate areas
	vector<double> areas;
	areas = area_triangles(points, triangles);

	for (auto _ : s) {
		sparse_mat B;
		B = assemble_matrix(points, triangles, lines, areas);
	}
}

//testing time evolution, exluding initial state and file writing, on sparse-man matrix B
// fixed mesh
// fixed time grid
static void timeEvolutionSparseCsrBench2(benchmark::State& s) {
	// Input data mesh
	vector<point> points;
	vector<triangle> triangles;
	vector<line> lines;
	read_mesh("assets/square.msh", points, triangles, lines);

	//Calculate areas
	vector<double> areas;
	areas = area_triangles(points, triangles);

	//assemble B
	sparse_mat B;
	B = assemble_matrix(points, triangles, lines, areas);

	//calculate initial state
	vector<double> ini_state(points.size());
	for (int i = 0; i < points.size(); i++) {
		double x = points[i].Getx(), y = points[i].Gety();
		ini_state[i] = heat(x, y);
	}

	for (auto _ : s) {
		vector<double> u = ini_state;
		int count = 0;
		double end_time = 0.1;
		double dt = 0.00001;
		int frame_count = int(end_time / dt);
		for (int i = 0; i < frame_count; i++) {
			u = one_timestep(dt, B, u); //calculates one timestep
			count += 1;
		}
	}
}


//testing how fast the sparse-man implentation is for a range of dt
//mesh is fixed
//initial heat condition is fixed
//begin and endtime are fixed
static void odeSparseCsrBench3(benchmark::State& s) {

	// Vector with different time steps dt to be benchmarked
	vector<double> steps(10);
	steps[0] = 0.00001;
	for (int i = 1; i < 10; i++) {
		steps[i] = steps[i - 1] + 0.00001;
	}

	// Input data mesh
	vector<point> points;
	vector<triangle> triangles;
	vector<line> lines;
	read_mesh("assets/square.msh", points, triangles, lines);

	//Calculate areas
	vector<double> areas;
	areas = area_triangles(points, triangles);

	for (auto _ : s) {
		s.PauseTiming();
		double dt = steps[s.range(0)];
		s.ResumeTiming();
		pdeSparseCsr(points, triangles, lines, areas, dt);
	}
}

// Register the benchmarks

BENCHMARK(assembleMatrixSparseCsrBench1)->Unit(benchmark::kMillisecond)->Iterations(30);
BENCHMARK(timeEvolutionSparseCsrBench2)->Unit(benchmark::kMillisecond)->Iterations(30);
BENCHMARK(odeSparseCsrBench3)->DenseRange(0, 9)->Unit(benchmark::kSecond)->Iterations(2);

// replaces normal main function
BENCHMARK_MAIN();
