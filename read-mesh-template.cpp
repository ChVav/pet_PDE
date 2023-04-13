#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

class point {
    double _x, _y;
public:
    point() {}
    point(double x, double y): _x(x), _y(y) {}
    double _Getx() {return _x;}
    double _Gety() {return _y;}
};

class line {
    int _i1, _i2;
public:
    line() {}
    line(int i1, int i2): _i1(i1), _i2(i2) {}
    int _GetPoint1() {return _i1;}
    int _GetPoint2() {return _i2;}
};

class triangle {
    int _i1, _i2, _i3;
public:
    triangle() {}
    triangle(int i1, int i2, int i3): _i1(i1), _i2(i2), _i3(i3) {}
    int _GetPoint1() {return _i1;}
    int _GetPoint2() {return _i2;}
    int _GetPoint3() {return _i3;}
    bool _has_vertex(int i)const {
        if (i==_i1 || i==_i2 || i==_i3){return true;}
        else {return false;}
    }
};

void read_mesh(string filename, vector<point>& points,
        vector<triangle>& triangles, vector<line>& lines) {
    ifstream fs(filename.c_str());
    string s="";

    // read until nodes
    // had to change from $Nodes to $ParametricNodes
    while( s != "$ParametricNodes")
        fs >> s;
    
    int num_nodes;
    fs >> num_nodes;
    points.resize(num_nodes);
    cout<<"number of nodes: "<<num_nodes<<endl;

    for(int i=0;i<num_nodes;i++) {
        int id;
        double x, y, unused;
        fs >> id >> x >> y >> unused;
        char c=' ';
        while( c != '\n' ) {
            fs.get(c);
        }
        // cout <<"i: "<<i<<", id "<<id<< ": x coord: "<<x<<"; y coord: "<<y<<endl;
        points[i]=point(x,y);
    }

    fs >> s;
    fs >> s;
    cout << s << " (should be $Elements)" << endl;
    //read in lines and triangles
    int num_elements;
    fs >> num_elements;
    //cout<< num_elements << endl;
    lines.reserve(num_elements);
    triangles.reserve(num_elements);
    /*
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
    }*/
    cout << endl;
}

int main(int argc, char* argv[])
{
    // read data
    vector<point> points;
    vector<triangle> triangles;
    vector<line> lines;
    read_mesh("square.msh", points, triangles, lines);

    cout << "Number of points:    " << points.size() << endl;
    //cout << "Number of triangles: " << triangles.size() << endl;
    //testing
    cout << "point 15 has coordinates x: "<<points[14]._Getx()<<" and y: "<<points[14]._Gety()<<endl;
    // TODO: write the .vtk file
}

