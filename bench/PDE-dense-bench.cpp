
/* Benchmarking dense implementation*/

#include <benchmark/benchmark.h>
#include "common.h"
#include "dense.h"

#include <cassert>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

//Functions unique to dense implementation

// calculate matrix B
vector<double> assemble_matrix(
	const vector<point>& points,
	const vector<triangle>& triangles,
	const vector<line>& lines,
	const vector<double>& areas) {

	int n_v = points.size();
	int n_T = triangles.size();

	// set B to zero
	vector<double> B(n_v * n_v);
	fill(begin(B), end(B), 0.0);

	// use pseudo code from finite_elements.pdf:
	for (int i = 0; i < n_v; i++) {
		double H = 0;
		for (int k = 0; k < n_T; k++) {
			if (triangles[k].has_vertex(i)) {
				H += areas[k];
			}
		}

		if (!on_boundary(i, lines)) {
			for (int j = 0; j < n_v; j++) {
				for (int k = 0; k < n_T; k++) {
					if (triangles[k].has_vertex(i) && triangles[k].has_vertex(j)) {
						double bik, bjk;
						double cik, cjk;
						compute_bc(i, triangles[k], points, bik, cik);
						compute_bc(j, triangles[k], points, bjk, cjk);
						B[j + n_v * i] -= areas[k] * (bik * bjk + cik * cjk);
					}
				}
				B[j + n_v * i] /= H;
			}
		}
	}
	return B;
}

// forward Euler, calculate one timestep
vector<double> one_timestep(double dt, const vector<double>& B, const vector<double>& u_n) {
	//returns u_n+1 = u_n+dt * B * u_n
	vector<double> u_n1(u_n.size());
	for (int i = 0; i < u_n.size(); i++) {
		u_n1[i] = u_n[i];
		for (int j = 0; j < u_n.size(); j++) {
			const auto index = j + i * u_n.size();
			const auto b = B[index];

			u_n1[i] += dt * b * u_n[j];

			assert(!std::isinf(u_n1[i]));
			assert(!std::isnan(u_n1[i]));
		}
	}
	return u_n1;
}

// Combined functions for benchmarking dense implementation (assemble_matrix + forward Euler + initial heat state): 
void pdeDense(const vector<point>& points,
	const vector<triangle>& triangles,
	const vector<line>& lines,
	const vector<double>& areas,
	double dt) {

	vector<double> B;
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

// testing assembly of dense matrix only
// fixed mesh
static void assembleMatrixDenseBench1(benchmark::State& s) {
	// Input data mesh
	vector<point> points;
	vector<triangle> triangles;
	vector<line> lines;
	read_mesh("assets/square.msh", points, triangles, lines);

	//Calculate areas
	vector<double> areas;
	areas = area_triangles(points, triangles);

	for (auto _ : s) {
		vector<double> B;
		B = assemble_matrix(points, triangles, lines, areas);
	}
}

//testing time evolution, exluding initial state and file writing, on dense matrix B
// fixed mesh
// fixed time grid
static void timeEvolutionDenseBench2(benchmark::State& s) {
	// Input data mesh
	vector<point> points;
	vector<triangle> triangles;
	vector<line> lines;
	read_mesh("assets/square.msh", points, triangles, lines);

	//Calculate areas
	vector<double> areas;
	areas = area_triangles(points, triangles);

	//assemble B
	vector<double> B;
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


//testing how fast the dense implentation is for a range of dt
//mesh is fixed
//initial heat condition is fixed
//begin and endtime are fixed
static void odeDenseBench3(benchmark::State& s) {

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
		pdeDense(points, triangles, lines, areas, dt);
	}
}

// Register the benchmarks

BENCHMARK(assembleMatrixDenseBench1)->Unit(benchmark::kMillisecond)->Iterations(30);
BENCHMARK(timeEvolutionDenseBench2)->Unit(benchmark::kMillisecond)->Iterations(30);
BENCHMARK(odeDenseBench3)->DenseRange(0, 9)->Unit(benchmark::kSecond)->Iterations(2);

// replaces normal main function
BENCHMARK_MAIN();
