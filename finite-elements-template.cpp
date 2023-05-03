#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string.h>


using namespace std;

class point {
private:
    double _x, _y, _z; //(x,y,z) coordinates of vertices in mesh with z=0
public:
    point(){}
    point(double x, double y, double z): _x(x), _y(y), _z(z) {}
    double _Getx() const {return _x;}
    double _Gety() const {return _y;}
    double _Getz() const {return _z; }
};

class line {
private:
    int _i1, _i2; //indices of vertices, two connecting points gives lines
public:
    line(){}
    line(int i1, int i2): _i1(i1), _i2(i2) {}
    int _GetPoint1() const {return _i1;}
    int _GetPoint2() const {return _i2;}
};

class triangle {
 private:
    int _i1, _i2, _i3; //indices of vertices, three connecting points gives triangles
public:
    triangle(){}
    triangle(int i1, int i2, int i3): _i1(i1), _i2(i2), _i3(i3) {}
    int _GetPoint1() const {return _i1;}
    int _GetPoint2() const {return _i2;}
    int _GetPoint3() const {return _i3;}
    bool _has_vertex(int i)const {
        if (i==_i1 || i==_i2 || i==_i3){return true;}
        else {return false;}
    } 
};

void read_mesh(string filename, vector<point>& points,
        vector<triangle>& triangles, vector<line>& lines) {
    ifstream fs(filename.c_str());
    string s="";

    /*
    if (!fs.good()) {
        __debugbreak();
    }
    */
    // read until nodes
    // had to change from $Nodes to $ParametricNodes
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


double area_triangle(point p1, point p2, point p3) {
    double area;
    double x1, y1, x2, y2, x3, y3;
    x1=p1._Getx();
    y1=p1._Gety();
    x2=p2._Getx();
    y2=p2._Gety();
    x3=p3._Getx();
    y3=p3._Gety();
    area=1.0/2*abs((x1-x3)*(y2-y1)-(x1-x2)*(y3-y1));
    return area;
    return 0;
}

vector<double> area_triangles(const vector<point>& points,
        const vector<triangle>& triangles) {
    //returns array with areas of each triangle
    vector<double> areas(triangles.size());
    for (int i=0; i<triangles.size(); i++){
        point p1, p2, p3;
        p1=points[triangles[i]._GetPoint1()];
        p2=points[triangles[i]._GetPoint2()];
        p3=points[triangles[i]._GetPoint3()];
        areas[i]=area_triangle(p1, p2, p3);
    }
    return areas;
}

void compute_bc(int i, triangle t, const vector<point>& points, double& b, double& c) {
    point p1,p2,p3;
    if(i == t._GetPoint1()) {
        p1 = points[t._GetPoint1()];
        p2 = points[t._GetPoint2()];
        p3 = points[t._GetPoint3()];
    } else if(i == t._GetPoint2()) {
        p1 = points[t._GetPoint2()];
        p2 = points[t._GetPoint1()];
        p3 = points[t._GetPoint3()];
    } else if(i == t._GetPoint3()) {
        p1 = points[t._GetPoint3()];
        p2 = points[t._GetPoint1()];
        p3 = points[t._GetPoint2()];
    } else {
        cout << "ERROR: vertex i is not part of triangle i" << endl;
        exit(1);
    }

    double xi=p1._Getx(),yi=p1._Gety();
    double xj=p2._Getx(),yj=p2._Gety();
    double xk=p3._Getx(),yk=p3._Gety();

    double norm=xi*yj-xi*yk-xj*yi+xj*yk+xk*yi-xk*yj;
    b=yj-yk;
    c=xk-xj;
    /* how to use compute_bc:
    double bik, bjk;
    double cik, cjk;
    compute_bc(i, triangles[k], points, bik, cik);
    compute_bc(j, triangles[k], points, bjk, cjk);
    */
}

bool on_boundary(int i, const vector<line>& lines) {
    for (int j=0; j<lines.size(); j++){
        if (lines[j]._GetPoint1()==i || lines[j]._GetPoint2()==i) {
            //cout << "point " << i << " is on boundary, see line " << j << endl;
            return true;
        } else {}
        
    }
    return false;
}

vector<double> assemble_matrix(const vector<point>& points,
        const vector<triangle>& triangles, const vector<line>& lines, const vector<double>& areas) {

    int n_v = points.size();
    int n_T = triangles.size();

    // set B to zero
    vector<double> B(n_v*n_v);
    fill(begin(B), end(B), 0.0);

    // use pseudo code from finite_elements.pdf:
    for (int i=0; i<n_v; i++) {
        double H=0;
        for (int k=0; k<n_T; k++) {
            if (triangles[k]._has_vertex(i)){
                H+= areas[k];
            }
        }

        if (!on_boundary(i, lines)){
            for (int j=0; j<n_v; j++) {
                for (int k=0; k<n_T; k++) {
                    if (triangles[k]._has_vertex(i) && triangles[k]._has_vertex(j)){
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
    /* how to use compute_bc:
    double bik, bjk;
    double cik, cjk;
    compute_bc(i, triangles[k], points, bik, cik);
    compute_bc(j, triangles[k], points, bjk, cjk);
    */
    return B;
}

double heat(double x, double y){
    // specify the function for initial condition
    double width=0.2;
    double res=1*exp(-pow(0.5-x, 2)/pow(width, 2) - pow(0.5-y,2)/pow(width, 2));
    //double res= pow(sin(x/width),2)+pow(cos(y/width),2);
    return res;
}

void output(double time_step, const vector<point>& points, const vector<line>& lines, const vector<triangle>& triangles, vector<double>& u){
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
        out << points[i]._Getx() << " " << points[i]._Gety() << " " << points[i]._Getz() << std::endl; // write values
    }
    out << "LINES " << lines.size() << " " << lines.size() * 3 << std::endl;
    for (int i = 0; i < lines.size(); i++) {
        out <<"2 "<< lines[i]._GetPoint1() << " " << lines[i]._GetPoint2() << std::endl; // write values
    }
    out << "TRIANGLE_STRIPS " << triangles.size() << " " << triangles.size() * 4 << std::endl;
    for (int i = 0; i < triangles.size(); i++) {
        out << "3 " << triangles[i]._GetPoint1() << " " << triangles[i]._GetPoint2() << " " << triangles[i]._GetPoint3() << std::endl; // write values
    }
    //header and geometry is written
    //now write the produced data:    
    out << "POINT_DATA "<< points.size() << endl;
    out << "SCALARS value double 1" << endl;
    out << "LOOKUP_TABLE default" <<endl;

    for (int i=0; i<points.size(); i++){
        double x=points[i]._Getx(), y=points[i]._Gety();
        out << u[i] << endl;
    }
    out.close();
}

vector<double> one_timestep(double dt, const vector<double>& B, vector<double>& u_n){
    //calculates one timestep
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

int main() {
    // read data
    vector<point> points;
    vector<triangle> triangles;
    vector<line> lines;
    read_mesh("square.msh", points, triangles, lines);

    cout << "Number of points:    " << points.size() << endl;
    cout << "Number of lines: " << lines.size() << endl;
    cout << "Number of triangles: " << triangles.size() << endl;

    vector<double> areas;
    areas=area_triangles(points, triangles);
    
    vector<double> B;
    B=assemble_matrix(points, triangles, lines, areas);
    
    
    // actual time evolution:
    //initial state
    vector<double> ini_state(points.size());
    for (int i=0; i<points.size(); i++){
        double x=points[i]._Getx(), y=points[i]._Gety();
        ini_state[i]=heat(x, y);
    }
    //write initial state into output file          
    output(0.0, points, lines, triangles, ini_state);
    
    vector<double> u=ini_state;
    int count=0;
    int file_count=1;
    double end_time=1000;
    double dt=0.1;
    int frame_count=int(end_time/dt);
    //saves 100 files in timeseries:
    int save_every=int(end_time/dt /100);
    for (int i=0; i<frame_count; i++) {
        //double dt=0.01;
        u=one_timestep(dt, B, u); //calculates one timestep
        count+=1;
        if (count==save_every) {
            output(file_count/100.0, points, lines, triangles, u);
            count=0;
            cout << "saved file " << file_count << endl;
            file_count+=1;
        }
    }
    
    /*
    //testing
    double total_area=0;
    for (int i=0; i<triangles.size(); i++) {
        total_area+=areas[i];
        //cout << "Area of triangle " << i << ": " << areas[i] << endl;
    }
    cout << "total area is: " << total_area << endl;
    int points_on_boundary=0;
    for (int i=0; i<points.size(); i++) {
        on_boundary(i, lines);
        if (on_boundary(i, lines)) {points_on_boundary +=1;};
    }
    cout << "points on boundary : "<< points_on_boundary << endl; 
    double b, c;
    compute_bc(460, triangles[0], points, b, c);
    cout << "b is: " << b << endl;
    cout << "c is: " << c << endl;

    int B_diff_f0_count=0;
    for (int i=0; i<B.size(); i++) {
        //cout << "B[i] is: " << B[i] << endl;
        if (B[i]!=0) {
            B_diff_f0_count+=1;
        }
    }
    cout << B_diff_f0_count << " -B matrix elements are different from zero" << endl;

    int count_all=0;
    int count_same=0;
    for (int i=0; i<points.size(); i++) {
        for (int j=0; j<points.size(); j++) {
            count_all+=1;
            if (B[j+i*points.size()]==B[j*points.size()+i]) {count_same+=1;}
        }
    }
    cout << "all elements: " << count_all << endl;
    cout << "same elements:" << count_same << endl;
    int diff_from_0_count=0;
    
    cout << "points greater than 0.1: " << diff_from_0_count << endl;
    */
    return 0;
}