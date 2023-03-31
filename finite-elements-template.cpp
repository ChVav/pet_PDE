#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;

class point {
    double _x, _y;
    // TODO: write functions and constructor(s)
};

class line {
    int _i1, _i2;
    // TODO: write functions and constructor(s)
};

class triangle {
    int _i1, _i2, _i3;
public:
    // TODO: write functions and constructor(s)
    bool has_vertex(int i)const {
        // TODO: return true if the triangle contains the vertex i
    } 
};

void read_mesh(string filename, vector<point>& points,
        vector<triangle>& triangles, vector<line>& lines) {
    // TODO
}

double area_triangle(point p1, point p2, point p3) {
    // TODO: return the area of the triangle with vertices p1, p2, and p3
}

vector<double> area_triangles(const vector<point>& points,
        const vector<triangle>& triangles) {

    vector<double> areas(triangles.size());

    // TODO: areas[i] should contain the area of triangles[i]

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
}

vector<double> assemble_matrix(const vector<point>& points,
        const vector<triangle>& triangles, const vector<line>& lines) {

    int n_v = points.size();
    int n_T = triangles.size();

    // set B to zero
    vector<double> B(n_v*n_v);
    fill(begin(B), end(B), 0.0);

    // TODO: assemble the matrix using the functions defined above.

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

