#include "common.h"

#include <cassert>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

// calculate bc
void compute_bc(
	int i,
	triangle t,
	const vector<point>& points,
	double& b,
	double& c) {
	point p1, p2, p3;
	if (i == t.GetPoint1()) {
		p1 = points[t.GetPoint1()];
		p2 = points[t.GetPoint2()];
		p3 = points[t.GetPoint3()];
	}
	else if (i == t.GetPoint2()) {
		p1 = points[t.GetPoint2()];
		p2 = points[t.GetPoint3()];
		p3 = points[t.GetPoint1()];
	}
	else if (i == t.GetPoint3()) {
		p1 = points[t.GetPoint3()];
		p2 = points[t.GetPoint1()];
		p3 = points[t.GetPoint2()];
	}
	else {
		cout << "ERROR: vertex i is not part of triangle i" << endl;
		exit(1);
	}

	double xi = p1.Getx(), yi = p1.Gety();
	double xj = p2.Getx(), yj = p2.Gety();
	double xk = p3.Getx(), yk = p3.Gety();

	
	//This is what was intended?
	//Kind of implodes, results in Inf and nan in the final .vtk files
	double norm = (xi * yj) - (xi * yk) - (xj * yi) + (xj * yk) + (xk * yi) - (xk * yj);
	//cout << "norm of " << i <<" ,"<<t.GetPoint2()<<", "<<t.GetPoint3()<<" is: "<< norm << endl;
	/*
	assert(!std::isnan(norm));
	assert(!std::isinf(norm));
	assert(abs(norm - 0.0) > std::numeric_limits<double>::epsilon());
	*/
	b = (yj - yk) / norm;
	c = (xk - xj) / norm;
	

	//b = (yj - yk)/(-0.001);
	//c = (xk - xj)/(-0.001);
}

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
	//returns u_n+1=u_n+dt*B*u_n
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

int main() {
	// read data
	vector<point> points;
	vector<triangle> triangles;
	vector<line> lines;
	vector<double> norms_of_bc;
	read_mesh("assets/square.msh", points, triangles, lines);

	cout << "Number of points:    " << points.size() << endl;
	cout << "Number of lines: " << lines.size() << endl;
	cout << "Number of triangles: " << triangles.size() << endl;

	vector<double> areas;
	areas = area_triangles(points, triangles);

	vector<double> B;
	B = assemble_matrix(points, triangles, lines, areas);

	// actual time evolution:
	//initial state
	vector<double> ini_state(points.size());
	for (int i = 0; i < points.size(); i++) {
		double x = points[i].Getx(), y = points[i].Gety();
		ini_state[i] = heat(x, y);
	}
	
	//write initial state into output file          
	output(0.0, points, lines, triangles, ini_state);

	//testing one timestep:
	/*
	vector<double> u = ini_state;
	u=one_timestep(0.01, B, ini_state);

	output(0.01, points, lines, triangles, u);
	*/
	
	vector<double> u = ini_state;
	int count = 0;
	int file_count = 1;
	double end_time = 0.1;
	double dt = 0.00001;
	int frame_count = int(end_time / dt);
	//saves 100 files in timeseries:
	int save_every = int(end_time / dt / 100);
	for (int i = 0; i < frame_count; i++) {
		//double dt=0.01;
		u = one_timestep(dt, B, u); //calculates one timestep
		count += 1;
		if (count == save_every) {
			output(file_count / 100.0, points, lines, triangles, u);
			count = 0;
			cout << "saved file " << file_count << endl;
			file_count += 1;
		}
	}
	
	return 0;
}