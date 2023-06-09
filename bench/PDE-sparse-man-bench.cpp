
/* Benchmarking dense implementation (assemble_matrix + forward Euler + initial heat state),
should probably test assemble_matrix and forward Euler separately*/

#include <benchmark/benchmark.h>
#include "sparse-man.h"

#include <cassert>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>

using namespace std;

// point definitions
point::point() {}
point::point(double x, double y, double z) :
	_x(x),
	_y(y),
	_z(z) {}

double point::Getx() const { return _x; }
double point::Gety() const { return _y; }
double point::Getz() const { return _z; }

// line definitions
line::line() {}
line::line(int i1, int i2) :
	_i1(i1),
	_i2(i2) {}

int line::GetPoint1() const { return _i1; }
int line::GetPoint2() const { return _i2; }

// triangle definitions
triangle::triangle() {}
triangle::triangle(int i1, int i2, int i3) :
	_i1(i1),
	_i2(i2),
	_i3(i3) {}

int triangle::GetPoint1() const { return _i1; }
int triangle::GetPoint2() const { return _i2; }
int triangle::GetPoint3() const { return _i3; }

bool triangle::has_vertex(int i)const {
	if (i == _i1 || i == _i2 || i == _i3) {
		return true;
	}
	else {
		return false;
	}
}

// definition function to read mesh
void read_mesh(string filename, vector<point>& points,
	vector<triangle>& triangles, vector<line>& lines) {
	ifstream fs(filename.c_str());
	string s = "";

	while (s != "$ParametricNodes") {
		fs >> s;
	}

	int num_nodes;
	fs >> num_nodes;
	points.resize(num_nodes);

	for (int i = 0;i < num_nodes;i++) {
		int id;
		double x, y, z, unused;
		fs >> id >> x >> y >> z >> unused;
		char c = ' ';
		while (c != '\n') {
			fs.get(c);
		}
		// cout <<"i: "<<i<<", id "<<id<< ": x coord: "<<x<<"; y coord: "<<y<< "; z coord: "<<z <<endl;
		points[i] = point(x, y, z);
	}

	fs >> s;
	fs >> s;

	//read in lines and triangles
	int num_elements;
	fs >> num_elements;
	//cout << num_elements << endl;
	lines.resize(num_elements);
	triangles.resize(num_elements);

	int line_count = 0, triangle_count = 0;
	for (int i = 0;i < num_elements;i++) {
		//cout << "in loop: "<<i<<endl;
		int id, type, unused1, unused2, unused3;
		fs >> id >> type >> unused1 >> unused2 >> unused3;
		// if the element is a line
		if (type == 1) {
			int i1, i2;
			fs >> i1 >> i2;
			lines[line_count] = line(i1 - 1, i2 - 1);
			line_count += 1;

			//if the element is a triangle
		}
		else if (type == 2) {
			int i1, i2, i3;
			fs >> i1 >> i2 >> i3;
			triangles[triangle_count] = triangle(i1 - 1, i2 - 1, i3 - 1);
			triangle_count += 1;
		}
		else {
			// read till the end of the line
			char c = ' ';
			while (c != '\n') {
				fs.get(c);
			}
		}
	}
	lines.resize(line_count);
	triangles.resize(triangle_count);
}

// definition function to calculate area of a triangle
double area_triangle(point p1, point p2, point p3) {
	double area;
	double x1, y1, x2, y2, x3, y3;
	x1 = p1.Getx();
	y1 = p1.Gety();
	x2 = p2.Getx();
	y2 = p2.Gety();
	x3 = p3.Getx();
	y3 = p3.Gety();
	area = 0.5 * abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
	return area;
}

// definition function to calculate areas of each triangle
vector<double> area_triangles(const vector<point>& points, const vector<triangle>& triangles) {
	//returns array with areas of each triangle
	vector<double> areas(triangles.size());
	for (int i = 0; i < triangles.size(); i++) {
		point p1, p2, p3;
		p1 = points[triangles[i].GetPoint1()];
		p2 = points[triangles[i].GetPoint2()];
		p3 = points[triangles[i].GetPoint3()];
		areas[i] = area_triangle(p1, p2, p3);
	}
	return areas;
}

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
		p2 = points[t.GetPoint1()];
		p3 = points[t.GetPoint3()];
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
	/*
	assert(!std::isnan(norm));
	assert(!std::isinf(norm));
	assert(abs(norm - 0.0) > std::numeric_limits<double>::epsilon());
	*/
	b = (yj - yk) / norm;
	c = (xk - xj) / norm;


	//b = (yj - yk);
	//c = (xk - xj);

}

// check boundary
bool on_boundary(int i, const vector<line>& lines) {
	for (int j = 0; j < lines.size(); j++) {
		if (lines[j].GetPoint1() == i || lines[j].GetPoint2() == i) {
			return true;
		}
		else {}
	}
	return false;
}

sparse_mat::sparse_mat() {}
sparse_mat::sparse_mat(int dim, vector<int> pos, vector<double> val) {
        _count=pos.size();
		_dim=dim;
		_pos=pos;
		_val=val;
    }


vector<double> sparse_mat::vector_mult(vector<double> u, double dt=1) {
    vector<double> res(u.size());
    for (int k=0; k<_count; k++) {
        const int i=_pos[k]/_dim;
        const int j=_pos[k]%_dim;
        res[i]+=dt*_val[k]*u[j];
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
				double Bij=0;
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
					positions.push_back(j+n_v*i);
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

// function for setting initial conditions
double heat(double x, double y) {
	double width = 0.2;
	double res = 1 * exp(-pow(0.5 - x, 2) / pow(width, 2) - pow(0.5 - y, 2) / pow(width, 2));
	//double res= pow(sin(x/width),2)+pow(cos(y/width),2);
	return res;
}

// one timestep
vector<double> one_timestep(double dt, sparse_mat& B, const vector<double>& u_n) {
	//returns u_n+1=u_n+dt*B*u_n
	vector<double> u_n1(u_n.size());
	vector<double> delta_u(u_n.size());
	delta_u=B.vector_mult(u_n, dt);
	for (int i=0; i<u_n.size(); i++) {
		u_n1[i]=delta_u[i]+u_n[i];
	}
	return u_n1;
}

// Final function to benchmark: PDE solver, dense implementation; not writing to output
void pdeDense(const vector<point>& points,
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
	//saves 100 files in timeseries:
	int save_every = int(end_time / dt / 100);
	for (int i = 0; i < frame_count; i++) {
		//double dt=0.01;
		u = one_timestep(dt, B, u); //calculates one timestep
		count += 1;
	}
}



// wrapper used by google benchmarking
//testing how fast the dense implentation is for a range of dt
//mesh is fixed
//initial heat condition is fixed
//begin and endtime are fixed
static void odeSparseManBench1(benchmark::State& s) {

	// Vector with different time steps dt to be benchmarked
	vector<double> steps(10);
	steps[0] = 0.000001;
	for (int i = 1; i < 10; i++) {
		steps[i] = steps[i - 1] + 0.000001;
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
		double dt = steps[s.range(0)];
		pdeDense(points, triangles, lines, areas, dt);
	}
}

// Register the benchmark, test a range of dt, usual input arguments are integers, here indices to the dt vector
BENCHMARK(odeSparseManBench1)->DenseRange(0, 9)->Unit(benchmark::kSecond)->Iterations(2);

// replaces normal main function
BENCHMARK_MAIN();
