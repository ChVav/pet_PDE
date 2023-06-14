
/* Benchmarking sparse-man implementation*/

#include <benchmark/benchmark.h>
#include "common.h"
#include "sparse-man.h"

#include <cassert>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>

using namespace std;

//Functions unique to sparse-man implementation
sparse_mat::sparse_mat() {};
sparse_mat::sparse_mat(int dim, vector<int> pos, vector<double> val) {
	_count = pos.size();
	_dim = dim;
	_pos = pos;
	_val = val;
}

vector<double> sparse_mat::vector_mult(vector<double> u, double dt = 1) {
	vector<double> res(u.size());
	for (int k = 0; k < _count; k++) {
		const int i = _pos[k] / _dim;
		const int j = _pos[k] % _dim;
		res[i] += dt * _val[k] * u[j];
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
	vector<int> positions;
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
					positions.push_back(j + n_v * i);
					values.push_back(Bij);
				}
			}
		}
	}
	// construct B
	sparse_mat B;
	B = sparse_mat(points.size(), positions, values);
	return B;
}

// forward Euler, calculate one timestep
vector<double> one_timestep(double dt, sparse_mat& B, const vector<double>& u_n) {
	//returns u_n+1 = u_n + dt * B * u_n
	vector<double> u_n1(u_n.size());
	vector<double> delta_u(u_n.size());
	delta_u = B.vector_mult(u_n, dt);
	for (int i = 0; i < u_n.size(); i++) {
		u_n1[i] = delta_u[i] + u_n[i];
	}
	return u_n1;
}

// Combined functions for benchmarking sparse-man implementation (assemble_matrix + forward Euler + initial heat state): 
void pdeSparseMan(const vector<point>& points,
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
static void assembleMatrixSparseManBench1(benchmark::State& s) {
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
static void timeEvolutionSparseManBench2(benchmark::State& s) {
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
static void odeSparseManBench3(benchmark::State& s) {

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
		pdeSparseMan(points, triangles, lines, areas, dt);
	}
}

// Register the benchmarks

BENCHMARK(assembleMatrixSparseManBench1)->Unit(benchmark::kMillisecond)->Iterations(30);
BENCHMARK(timeEvolutionSparseManBench2)->Unit(benchmark::kMillisecond)->Iterations(30);
BENCHMARK(odeSparseManBench3)->DenseRange(0, 9)->Unit(benchmark::kSecond)->Iterations(2);

// replaces normal main function
BENCHMARK_MAIN();
