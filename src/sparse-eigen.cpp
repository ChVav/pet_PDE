#include "common.h"
#include "sparse-eigen.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <cassert>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

// calculate matrix B
Eigen::SparseMatrix<double> assemble_matrix(
	const vector<point>& points,
	const vector<triangle>& triangles,
	const vector<line>& lines,
	const vector<double>& areas) {

	int n_v = points.size();
	int n_T = triangles.size();

	Eigen::SparseMatrix<double> B(n_v, n_v);

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
					B.insert(i, j) = Bij;
					B.makeCompressed();
				}
			}
		}
	}
	return B;
}

// forward Euler, calculate one timestep, for first step, making a copy to standard vector format (should rewrite more off the now common functions to do this more efficiently)
Eigen::VectorXd one_timestep(double dt, Eigen::SparseMatrix<double>& B, Eigen::VectorXd u_n) {
	auto u_n1 = u_n + dt * (B * u_n);
	return u_n1;
}

//Overloaded output function
void output(
	double time_step,
	const vector<point>& points,
	const vector<line>& lines,
	const vector<triangle>& triangles,
	Eigen::VectorXd u) {

	// write the .vtk file, output is vector<double> u
	string file1 = "mesh_t_";
	string file2 = to_string(time_step);
	string file3 = ".vtk";
	string filename = "data/" + file1 + file2 + file3;
	ofstream out(filename); // format refer to https://kitware.github.io/vtk-examples/site/VTKFileFormats/
	out << "# vtk DataFile Version 3.0" << endl;// Header: file version and identifier
	out << "Mesh" << endl;// Title
	out << "ASCII" << endl;// Data type (file format)
	out << "DATASET POLYDATA" << endl;// dataset structure, think should be polytonal data
	out << "POINTS " << points.size() << " double" << endl;
	for (int i = 0; i < points.size(); i++) {
		out << points[i].Getx() << " " << points[i].Gety() << " " << points[i].Getz() << endl; // write values
	}
	out << "LINES " << lines.size() << " " << lines.size() * 3 << endl;
	for (int i = 0; i < lines.size(); i++) {
		out << "2 " << lines[i].GetPoint1() << " " << lines[i].GetPoint2() << endl; // write values
	}
	out << "TRIANGLE_STRIPS " << triangles.size() << " " << triangles.size() * 4 << endl;
	for (int i = 0; i < triangles.size(); i++) {
		out << "3 " << triangles[i].GetPoint1() << " " << triangles[i].GetPoint2() << " " << triangles[i].GetPoint3() << endl; // write values
	}
	//header and geometry is written
	//now write the produced data:
	out << "POINT_DATA " << points.size() << endl;
	out << "SCALARS value double 1" << endl;
	out << "LOOKUP_TABLE default" << endl;

	for (int i = 0; i < points.size(); i++) {
		double x = points[i].Getx(), y = points[i].Gety();
		out << u[i] << endl;
	}
	out.close();
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

	Eigen::SparseMatrix<double> B = assemble_matrix(points, triangles, lines, areas);

	// actual time evolution:
	//initial state
	vector<double> ini_state(points.size());
	for (int i = 0; i < points.size(); i++) {
		double x = points[i].Getx(), y = points[i].Gety();
		ini_state[i] = heat(x, y);
	}

	//create a pointer for initial state to share memory with an Eigen::Vector
	Eigen::VectorXd u = Eigen::Map<Eigen::VectorXd>(ini_state.data(), ini_state.size());

	//write initial state into output file
	output(0.0, points, lines, triangles, u);

	//time evolution
	int count = 0;
	int file_count = 1;
	double end_time = 0.1;
	double dt = 0.00001;
	int frame_count = int(end_time / dt);
	//saves 100 files in timeseries:
	int save_every = int(end_time / dt / 100);

	for (int i = 0; i < frame_count; i++) {
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