#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;

class point {
    double _x, _y;
    // TODO: write functions and constructor(s)
public:
    point() : _x(0.0), _y(0.0) {}
    point(double x, double y) : _x(x), _y(y) {}
    double x() const { return _x; }
    double y() const { return _y; }
};

class line {
    int _i1, _i2;
    // TODO: write functions and constructor(s)
public:
    line() : _i1(-1), _i2(-1) {}
    line(int i1, int i2) : _i1(i1), _i2(i2) {}
    int i1() const { return _i1; }
    int i2() const { return _i2; }
};

class triangle {
    int _i1, _i2, _i3;
public:
    // TODO: write functions and constructor(s)
    triangle() : i1(0), i2(0), i3(0) {}
    triangle(int i1, int i2, int i3) : _i1(i1), _i2(i2), _i3(i3) {}
    int i1() const { return _i1; }
    int i2() const { return _i2; }
    int i3() const { return _i3; }

    bool has_vertex(int i)const {
        // TODO: return true if the triangle contains the vertex i
        if (i==_i1 || i==_i2 || i==_i3){return true;}
        else {return false;}
    } 
};

void read_mesh(string filename, vector<point>& points,
        vector<triangle>& triangles, vector<line>& lines) {
    // TODO
    if (!file) {
    cerr << "Error opening file " << filename << endl;
    exit(1);
    }

    string line;
    int n_points, n_triangles, n_lines;
    getline(file, line);
    file >> n_points >> n_triangles >> n_lines;

    points.resize(n_points);
    for (int i = 0; i < n_points; ++i) {
        double x, y;
        file >> x >> y;
        points[i] = point(x, y);
    }

    triangles.resize(n_triangles);
    for (int i = 0; i < n_triangles; ++i) {
        int i1, i2, i3;
        file >> i1 >> i2 >> i3;
        triangles[i] = triangle(i1, i2, i3);
    }

    lines.resize(n_lines);
    for (int i = 0; i < n_lines; ++i) {
        int i1, i2;
        file >> i1 >> i2;
        lines[i] = line(i1, i2);
    }
}

double area_triangle(point p1, point p2, point p3) {
    // TODO: return the area of the triangle with vertices p1, p2, and p3
    double area = 0.5 * fabs((p2.x - p1.x) * (p3.y - p1.y) - (p3.x - p1.x) * (p2.y - p1.y));
    return area;

}

vector<double> area_triangles(const vector<point>& points,
        const vector<triangle>& triangles) {

    vector<double> areas(triangles.size());

    // TODO: areas[i] should contain the area of triangles[i]
    for (int i = 0; i < triangles.size(); i++) {
        int i1 = triangles[i].i1;
        int i2 = triangles[i].i2;
        int i3 = triangles[i].i3;
        point p1 = points[i1];
        point p2 = points[i2];
        point p3 = points[i3];
        areas[i] = area_triangle(p1, p2, p3);
    }


    return areas;
}

void compute_bc(int i, triangle t, const vector<point>& points, double& b, double& c) {
    point p1,p2,p3;
    if(i == t.i1) {
        p1 = points[t.i1];
        p2 = points[t.i2];
        p3 = points[t.i3];
    } else if(i == t.i2) {
        p1 = points[t.i2];
        p2 = points[t.i1];
        p3 = points[t.i3];
    } else if(i == t.i3) {
        p1 = points[t.i3];
        p2 = points[t.i1];
        p3 = points[t.i2];
    } else {
        cout << "ERROR: vertex i is not part of triangle t" << endl;
        exit(1);
    }

    // TODO: compute the coefficients b and c (i.e. b_{ik} and c_{ik}) from the lecture.
    // The code above makes sure that p1 always corresponds to vertex i.
}

bool on_boundary(int i, const vector<line>& lines) {

    // TODO: retun true if the vertex with index i is a boundary point, false // otherwise.
    for (int j = 0; j < lines.size(); j++) {
        if (lines[j]._i1 == i || lines[j]._i2 == i) {
            return true;
        }
    }
    return false;
}

vector<double> assemble_matrix(const vector<point>& points,
        const vector<triangle>& triangles, const vector<line>& lines) {

    int n_v = points.size();
    int n_T = triangles.size();

    // set B to zero
    vector<double> B(n_v*n_v);
    fill(begin(B), end(B), 0.0);

    // TODO: assemble the matrix using the functions defined above.
    for (int i = 0; i < n_v; i++) {
        for (int j = 0; j < n_v; j++) {
            if (i == j) {
                // diagonal entry
                double area_sum = 0.0;
                for (int k = 0; k < n_T; k++) {
                    if (triangles[k].has_vertex(i)) {
                        double area_k = area_triangle(points[triangles[k].i1],
                                points[triangles[k].i2], points[triangles[k].i3]);
                        area_sum += area_k / 3.0;
                    }
                }
                B[i*n_v + j] = area_sum;
            } else {
                // off-diagonal entry
                double area_sum = 0.0;
                for (int k = 0; k < n_T; k++) {
                    if (triangles[k].has_vertex(i) && triangles[k].has_vertex(j)) {
                        double b_ij, c_ij;
                        compute_bc(i, triangles[k], points, b_ij, c_ij);
                        double area_k = area_triangle(points[triangles[k].i1],
                                points[triangles[k].i2], points[triangles[k].i3]);
                        area_sum += b_ij * c_ij * area_k / 12.0;
                    }
                }
                B[i*n_v + j] = area_sum;
            }
        }
    }


    return B;
}

int main() {
    // read data
    vector<point> points;
    vector<triangle> triangles;
    vector<line> lines;
    read_mesh("test.msh", points, triangles, lines);

    cout << "Number of points:    " << points.size() << endl;
    cout << "Number of triangles: " << triangles.size() << endl;


    // TODO: assemble the matrix and solve the equation
}

