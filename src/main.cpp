#include "dense.h"

#include <iostream>

using namespace std;

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
        double x=points[i].Getx(), y=points[i].Gety();
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
    
    return 0;
}