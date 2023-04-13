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
    //cout << num_elements << endl;
    lines.reserve(num_elements);
    triangles.reserve(num_elements);
    
    cout << "read element of type "<<endl;
    
    int line_count=0, triangle_count=0;
    for(int i=0;i<num_elements;i++) {
        //cout << "in loop: "<<i<<endl;
        int id, type, unused1, unused2, unused3;
        fs >> id >> type >> unused1 >> unused2 >> unused3;
        // if the element is a line
        if (type==1){
            int i1, i2;
            fs >> i1 >> i2;
            lines[line_count]=line(i1,i2);
            line_count+=1;
            //cout << "line added" << endl;
            //cout << i1 << " " << i2 << endl;
        
        //if the element is a triangle
        } else if (type==2) {
            int i1, i2, i3;
            fs >> i1 >> i2 >> i3;
            triangles[triangle_count]=triangle(i1,i2,i3);
            triangle_count+=1;
            //cout << "triangle added"<<endl;
            //cout << i1 <<" "<<i2<<" "<<i3<< endl;
        } else {
            // read till the end of the line
            //cout << "nothing added" << endl;
            //cout << id <<" "<< type <<" "<<unused1 <<" "<<unused2 <<" "<<unused3<< endl;
            char c=' ';
            while( c != '\n' ) {
                fs.get(c);
            }
        }
    }
    lines.resize(line_count);
    triangles.resize(triangle_count);
}

int main(int argc, char* argv[])
{
    // read data
    vector<point> points;
    vector<triangle> triangles;
    vector<line> lines;
    read_mesh("square.msh", points, triangles, lines);

    cout << "Number of points:    " << points.size() << endl;
    cout << "Number of lines: " << lines.size() << endl;
    cout << "Number of triangles: " << triangles.size() << endl;
    //testing
    cout << "point 15 has coordinates x: "<<points[14]._Getx()<<" and y: "<<points[14]._Gety()<<endl;
    cout << "coordinates of point 5: " << points[4]._Getx()<<" "<<points[4]._Gety()<<endl;
    cout << "coordinates of point 6: " << points[5]._Getx()<<" "<<points[5]._Gety()<<endl;
    cout << "coordinates of point 7: " << points[6]._Getx()<<" "<<points[6]._Gety()<<endl;
    cout << "line 0 has the points: " << lines[0]._GetPoint1()<<" "<<lines[0]._GetPoint2()<<endl;
    cout << "line 1 has the points: " << lines[1]._GetPoint1()<<" "<<lines[1]._GetPoint2()<<endl;
    cout << "triangle 0 has the points: "<<triangles[0]._GetPoint1()<<" "<<triangles[0]._GetPoint2()<<" "<<triangles[0]._GetPoint3()<<endl;
    cout << "triangle 0 contains point 3: " <<triangles[0]._has_vertex(3)<< endl;
    cout << "triangle 0 contains point 461: " <<triangles[0]._has_vertex(461)<< endl;
    cout << "triangle 0 contains point 391: " <<triangles[0]._has_vertex(391)<< endl;
    cout << "triangle 0 contains point 493: " <<triangles[0]._has_vertex(493)<< endl;
   
    
    // TODO: write the .vtk file
}

