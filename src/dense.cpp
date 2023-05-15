#include "dense.h"

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

// point definitions
point::point(){}
point::point(double x, double y, double z): 
    _x(x), 
    _y(y), 
    _z(z) {}

double point::Getx() const {return _x;}
double point::Gety() const {return _y;}
double point::Getz() const {return _z;}

// line definitions
line::line(){}
line::line(int i1, int i2): 
    _i1(i1), 
    _i2(i2) {}

int line::GetPoint1() const {return _i1;}
int line::GetPoint2() const {return _i2;}

// triangle definitions
triangle::triangle(){}
triangle::triangle(int i1, int i2, int i3): 
    _i1(i1), 
    _i2(i2), 
    _i3(i3) {}

int triangle::GetPoint1() const {return _i1;}
int triangle::GetPoint2() const {return _i2;}
int triangle::GetPoint3() const {return _i3;}

bool triangle::has_vertex(int i)const {
        if (i==_i1 || i==_i2 || i==_i3){
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
    string s="";

    while (s != "$ParametricNodes") {
        fs >> s;
    }
    
    int num_nodes;
    fs >> num_nodes;
    points.resize(num_nodes);
    cout<<"number of nodes: "<<num_nodes<<endl;

    for(int i=0;i<num_nodes;i++) {
        int id;
        double x, y, z, unused;
        fs >> id >> x >> y >> z >> unused;
        char c=' ';
        while( c != '\n' ) {
            fs.get(c);
        }
        // cout <<"i: "<<i<<", id "<<id<< ": x coord: "<<x<<"; y coord: "<<y<< "; z coord: "<<z <<endl;
        points[i]=point(x,y,z);
    }

    fs >> s;
    fs >> s;
    cout << s << " (should be $Elements)" << endl;
    //read in lines and triangles
    int num_elements;
    fs >> num_elements;
    //cout << num_elements << endl;
    lines.resize(num_elements);
    triangles.resize(num_elements);
    
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
            lines[line_count]=line(i1-1,i2-1);
            line_count+=1;
            //cout << "line added" << endl;
            //cout << i1 << " " << i2 << endl;
        
        //if the element is a triangle
        } else if (type==2) {
            int i1, i2, i3;
            fs >> i1 >> i2 >> i3;
            triangles[triangle_count]=triangle(i1-1,i2-1,i3-1);
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

// definition function to calculate area of a triangle
double area_triangle(point p1, point p2, point p3) {
    double area;
    double x1, y1, x2, y2, x3, y3;
    x1=p1.Getx();
    y1=p1.Gety();
    x2=p2.Getx();
    y2=p2.Gety();
    x3=p3.Getx();
    y3=p3.Gety();
    area=1.0/2*abs((x1-x3)*(y2-y1)-(x1-x2)*(y3-y1));
    return area;
}

// definition function to calculate areas of each triangle
vector<double> area_triangles(const vector<point>& points, const vector<triangle>& triangles) {
    //returns array with areas of each triangle
    vector<double> areas(triangles.size());
    for (int i=0; i<triangles.size(); i++){
        point p1, p2, p3;
        p1=points[triangles[i].GetPoint1()];
        p2=points[triangles[i].GetPoint2()];
        p3=points[triangles[i].GetPoint3()];
        areas[i]=area_triangle(p1, p2, p3);
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
    point p1,p2,p3;
    if(i == t.GetPoint1()) {
        p1 = points[t.GetPoint1()];
        p2 = points[t.GetPoint2()];
        p3 = points[t.GetPoint3()];
    } else if(i == t.GetPoint2()) {
        p1 = points[t.GetPoint2()];
        p2 = points[t.GetPoint1()];
        p3 = points[t.GetPoint3()];
    } else if(i == t.GetPoint3()) {
        p1 = points[t.GetPoint3()];
        p2 = points[t.GetPoint1()];
        p3 = points[t.GetPoint2()];
    } else {
        cout << "ERROR: vertex i is not part of triangle i" << endl;
        exit(1);
    }

    double xi=p1.Getx(),yi=p1.Gety();
    double xj=p2.Getx(),yj=p2.Gety();
    double xk=p3.Getx(),yk=p3.Gety();

    double norm=xi*yj-xi*yk-xj*yi+xj*yk+xk*yi-xk*yj;
    b=yj-yk;
    c=xk-xj;
}

// check boundary
bool on_boundary(int i, const vector<line>& lines) {
    for (int j=0; j<lines.size(); j++){
        if (lines[j].GetPoint1()==i || lines[j].GetPoint2()==i) {
            //cout << "point " << i << " is on boundary, see line " << j << endl;
            return true;
        } else {}
        
    }
    return false;
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
    vector<double> B(n_v*n_v);
    fill(begin(B), end(B), 0.0);

    // use pseudo code from finite_elements.pdf:
    for (int i=0; i<n_v; i++) {
        double H=0;
        for (int k=0; k<n_T; k++) {
            if (triangles[k].has_vertex(i)){
                H+= areas[k];
            }
        }

        if (!on_boundary(i, lines)){
            for (int j=0; j<n_v; j++) {
                for (int k=0; k<n_T; k++) {
                    if (triangles[k].has_vertex(i) && triangles[k].has_vertex(j)){
                        double bik, bjk;
                        double cik, cjk;
                        compute_bc(i, triangles[k], points, bik, cik);
                        compute_bc(j, triangles[k], points, bjk, cjk);
                        B[j+ n_v*i] -= areas[k]*(bik*bjk+cik*cjk);
                    }
                }
                B[j+n_v*i] /= H;
            }
        }
    }
    return B;
}

// function for setting initial conditions
double heat(double x, double y){
    double width=0.2;
    double res=1*exp(-pow(0.5-x, 2)/pow(width, 2) - pow(0.5-y,2)/pow(width, 2));
    //double res= pow(sin(x/width),2)+pow(cos(y/width),2);
    return res;
}

// function to write .vtk file
void output(
    double time_step,
    const vector<point>& points,
    const vector<line>& lines,
    const vector<triangle>& triangles,
    vector<double>& u) {

    // write the .vtk file, output is vector<double> u
    string file1 = "mesh_t_";
    string file2 = to_string(time_step);
    string file3 =".vtk";
    string filename= "data/" + file1 + file2 + file3;
    std::ofstream out(filename); // format refer to https://kitware.github.io/vtk-examples/site/VTKFileFormats/
    out << "# vtk DataFile Version 3.0" << std::endl;// Header: file version and identifier
    out << "Mesh" << std::endl;// Title
    out << "ASCII" << std::endl;// Data type (file format)
    out << "DATASET POLYDATA" << std::endl;// dataset structure, think should be polytonal data
    out << "POINTS " << points.size() <<" double" << std::endl;
    for (int i = 0; i < points.size(); i++) {
        out << points[i].Getx() << " " << points[i].Gety() << " " << points[i].Getz() << std::endl; // write values
    }
    out << "LINES " << lines.size() << " " << lines.size() * 3 << std::endl;
    for (int i = 0; i < lines.size(); i++) {
        out <<"2 "<< lines[i].GetPoint1() << " " << lines[i].GetPoint2() << std::endl; // write values
    }
    out << "TRIANGLE_STRIPS " << triangles.size() << " " << triangles.size() * 4 << std::endl;
    for (int i = 0; i < triangles.size(); i++) {
        out << "3 " << triangles[i].GetPoint1() << " " << triangles[i].GetPoint2() << " " << triangles[i].GetPoint3() << std::endl; // write values
    }
    //header and geometry is written
    //now write the produced data:    
    out << "POINT_DATA "<< points.size() << endl;
    out << "SCALARS value double 1" << endl;
    out << "LOOKUP_TABLE default" <<endl;

    for (int i=0; i<points.size(); i++){
        double x=points[i].Getx(), y=points[i].Gety();
        out << u[i] << endl;
    }
    out.close();
}

// forward Euler, calculate one timestep
vector<double> one_timestep(double dt, const vector<double>& B, const vector<double>& u_n){
    //returns u_n+1=u_n+dt*B*u_n
    vector<double> u_n1(u_n.size());
    for (int i=0; i<u_n.size(); i++){
        u_n1[i]=u_n[i];
        for (int j=0; j<u_n.size(); j++) {
            u_n1[i]+=dt*B[j+i*u_n.size()]*u_n[j];
        }
    }
    return u_n1;
}