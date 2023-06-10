#pragma once

#include <vector>
#include <string>

class point
{ // stores (x,y,z) coordinates of vertices in mesh with z=0
public:
	// getters
	double Getx() const;
	double Gety() const;
	double Getz() const;

	// constructor, can we get rid of default constructors??
	point();
	point(double x, double y, double z);

private:
	double _x, _y, _z;
};

class line
{ //stores indices of vertices, two connecting points gives lines
public:
	int GetPoint1() const;
	int GetPoint2() const;

	line();
	line(int i1, int i2);

private:
	int _i1, _i2;
};

class triangle
{ //stores indices of vertices, three connecting points gives triangles
public:
	bool has_vertex(int i)const;

	int GetPoint1() const;
	int GetPoint2() const;
	int GetPoint3() const;

	triangle();
	triangle(int i1, int i2, int i3);

private:
	int _i1, _i2, _i3;
};

void read_mesh(std::string filename, std::vector<point>& points, std::vector<triangle>& triangles, std::vector<line>& lines);

double area_triangle(point p1, point p2, point p3);

std::vector<double> area_triangles(const std::vector<point>& points, const std::vector<triangle>& triangles);

void compute_bc(int i, triangle t, const std::vector<point>& points, double& b, double& c);

bool on_boundary(int i, const std::vector<line>& lines);

double heat(double x, double y);

void output(
	double time_step,
	const std::vector<point>& points,
	const std::vector<line>& lines,
	const std::vector<triangle>& triangles,
	std::vector<double>& u);

std::vector<double> one_timestep(
	double dt,
	const std::vector<double>& B,
	const std::vector<double>& u_n);
