#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

class point {
private:
    double _x, _y;
    // TODO: write functions and constructor(s)
};

class line {
private:
    int _i1, _i2;
    // TODO: write functions and constructor(s)
};

class triangle {
private:
    int _i1, _i2, _i3;
    // TODO: write functions and constructor(s)
};

void read_mesh(string filename, vector<point>& points,
        vector<triangle>& triangles, vector<line>& lines) {
    ifstream fs(filename.c_str());
    string s="";

    // read until nodes
    while( s != "$Nodes")
        fs >> s;

    int num_nodes;
    fs >> num_nodes;

    for(int i=0;i<num_nodes;i++) {
        int id;
        double x, y, unused;
        fs >> id >> x >> y >> unused;

        // TODO: read data and populate the vector points
    }

    fs >> s;
    fs >> s;
    cout << s << " (should be $Elements)" << endl;

    int num_elements;
    fs >> num_elements;

    cout << "read element of type ";
    for(int i=0;i<num_elements;i++) {
        int id, type, unused;
        fs >> id >> type >> unused >> unused >> unused;
        cout << type << " ";

        // TODO: check the type (either 1 for line or 2 for triangle) and
        // populate the vectors lines and triangles.
        } else {
            // read till the end of the line
            char c=' ';
            while( c != '\n' ) {
                fs.get(c);
            }
        }
    }
    cout << endl;
}

int main(int argc, char* argv[])
{
    // read data
    vector<point> points;
    vector<triangle> triangles;
    vector<line> lines;
    read_mesh("test.msh", points, triangles, lines);

    cout << "Number of points:    " << points.size() << endl;
    cout << "Number of triangles: " << triangles.size() << endl;

    // TODO: write the .vtk file
}

