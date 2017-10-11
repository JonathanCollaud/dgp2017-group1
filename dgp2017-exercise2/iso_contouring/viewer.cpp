//=============================================================================
//
//   Code framework for the lecture
//
//   "Digital 3D Geometry Processing"
//
//   Copyright (C) 2017 by Computer Graphics and Geometry Laboratory,
//         EPF Lausanne
//
//-----------------------------------------------------------------------------
#include "viewer.h"
#include <surface_mesh/Surface_mesh.h>

using std::min;
using std::max;
using namespace surface_mesh;


using std::cout ;using std::endl ;


typedef Surface_mesh Mesh;

// ========================================================================
// NOTE : We've only included the functions you need to implement (or the
//        ones that you will need to use) in the cpp file. This is not the
//        best practice as you normaly would have all the implementation of
//        the functions here and only the declaration in the header file
//        but it allows you to have all the things you need here.
// ========================================================================

// ========================================================================
// EXERCISE 1
// ========================================================================
Scalar Viewer::iso_value(Point v_pos)
{
    float x,y;
    x = v_pos.x;
    y = v_pos.y;

    // ----- (un)comment a line to change the function you are valing
    //Scalar iso = sqrt(x*x + y*y) - 1;
    //Scalar iso = sin(2*x+2*y) - cos(4*x*y) +1;
    Scalar iso = y*y - sin(x*x);

    return iso;
}

void Viewer::calc_iso_contouring() {
    Mesh::Vertex_property<Scalar> v_iso = mesh.vertex_property<Scalar>("v:iso", 0);
    segment_points.clear();
    std::vector<Point> v_positions(mesh.n_vertices());
    std::vector<std::vector<int> > triangle_ids;

    for (auto v: mesh.vertices()) {
        Point v_pos = mesh.position(v);
        v_positions[v.idx()] = v_pos;
        Scalar iso = 0;

        iso = iso_value(v_pos);

        v_iso[v] = iso; //this variable is for coloring the density; do not change this variable
    }

    for(auto f: mesh.faces()) {
        std::vector<int> vv(3);
        int k = 0;
        for (auto v: mesh.vertices(f)) {
            //cout<<mesh.vertices(f)<< endl ;
            vv[k] = v.idx();
            ++k;
        }
        triangle_ids.push_back(vv);
    }

    //segment_points is defined in viewer.h as std::vector<Point> segment_points;
    //add points in segment_points forming an edge one after the other;
    //for example segment_points[0] and segment_points[1] are two points forming the first edge
    //and segment_points[2] and segment_points[3] are two points forming the second edge

    // ----- add your code here -----

    for(auto f: mesh.faces()) {
        Vec3 val, x, y, x0 = {0,0,0}, y0 = {0,0,0};
        int i = 0;
        int intersection_count = 0;

        // Compute value on each vertex of the face
        for (auto v: mesh.vertices(f)) {
            x[i] = v_positions[v.idx()][0] ;
            y[i] = v_positions[v.idx()][1] ;
            val[i] = iso_value(v_positions[v.idx()]) ;
            i++ ;
        }

        // Found intersections
        for(i = 0; i < 3; i++){
            if(val[i] == 0){
                x0[intersection_count] = x[i];
                y0[intersection_count] = y[i];
                intersection_count++;
            }
            if(val[i]*val[(i+1)%3] < 0){
                x0[intersection_count] = x[i] + (val[i]/(val[i]-val[(i+1)%3]))*(x[(i+1)%3]-x[i]);
                y0[intersection_count] = y[i] + (val[i]/(val[i]-val[(i+1)%3]))*(y[(i+1)%3]-y[i]);
                intersection_count++;
            }
        }

        // Add intersections to segment_points
        for(i = 0; i <=3; i++){
            segment_points.push_back({x0[i],y0[i],0});
        }
    }

    // ------------------------------
}
