//=============================================================================
//
//   Code framework for the lecture
//
//   "Digital 3D Geometry Processing"
//
//   Gaspard Zoss
//
//   Copyright (C) 2016 by Computer Graphics and Geometry Laboratory,
//         EPF Lausanne
//
//   Edited 2017
//-----------------------------------------------------------------------------
#include "viewer.h"
#include <surface_mesh/Surface_mesh.h>

using std::min;
using std::max;
using namespace surface_mesh;

typedef Surface_mesh Mesh;

// ========================================================================
// NOTE : We've only included the functions you need to implement (or the
//        ones that you will need to use) in the cpp file. This is not the
//        best practice as you normaly would have all the implementation of
//        the functions here and only the declaration in the header file
//        but it allows you to have all the things you need here.
// ========================================================================

void Viewer::calc_weights() {
    calc_edges_weights();
    calc_vertices_weights();
}

void Viewer::calc_edges_weights() {
    Mesh::Halfedge h0, h1, h2;
    Mesh::Vertex v0, v1;
    Point p0, p1, p2, d0, d1;
    Scalar w;
    auto eweight = mesh.edge_property<Scalar>("e:weight", 0);
    for (auto e: mesh.edges()) {
        w = 0.0;

        h0 = mesh.halfedge(e, 0);
        v0 = mesh.to_vertex(h0);
        p0 = mesh.position(v0);

        h1 = mesh.halfedge(e, 1);
        v1 = mesh.to_vertex(h1);
        p1 = mesh.position(v1);

        h2 = mesh.next_halfedge(h0);
        p2 = mesh.position(mesh.to_vertex(h2));
        d0 = normalize(p0 - p2);
        d1 = normalize(p1 - p2);
        w += 1.0 / tan(acos(min(0.99f, max(-0.99f, dot(d0, d1)))));

        h2 = mesh.next_halfedge(h1);
        p2 = mesh.position(mesh.to_vertex(h2));
        d0 = normalize(p0 - p2);
        d1 = normalize(p1 - p2);
        w += 1.0 / tan(acos(min(0.99f, max(-0.99f, dot(d0, d1)))));

        w = max(0.0f, w);
        eweight[e] = w * 0.5;
    }
}

void Viewer::calc_vertices_weights() {
    Mesh::Face_around_vertex_circulator vf_c, vf_end;
    Mesh::Vertex_around_face_circulator fv_c;
    Scalar area;
    auto vweight = mesh.vertex_property<Scalar>("v:weight", 0);

    for (auto v: mesh.vertices()) {
        area = 0.0;
        vf_c = mesh.faces(v);

        if(!vf_c) {
            continue;
        }

        vf_end = vf_c;

        do {
            fv_c = mesh.vertices(*vf_c);

            const Point& P = mesh.position(*fv_c);  ++fv_c;
            const Point& Q = mesh.position(*fv_c);  ++fv_c;
            const Point& R = mesh.position(*fv_c);

            area += norm(cross(Q-P, R-P)) * 0.5f * 0.3333f;

        } while(++vf_c != vf_end);

        vweight[v] = 0.5 / area;
    }
}

void Viewer::computeValence() {
    Mesh::Vertex_property<Scalar> vertex_valence =
            mesh.vertex_property<Scalar>("v:valence", 0);
    for (auto v: mesh.vertices()) {
        vertex_valence[v] = mesh.valence(v);
    }
}
// ========================================================================
// EXERCISE 1.1
// ========================================================================
void Viewer::computeNormalsWithConstantWeights() {
    Point default_normal(0.0, 1.0, 0.0);
    Mesh::Vertex_property<Point> v_cste_weights_n =
            mesh.vertex_property<Point>("v:cste_weights_n", default_normal);

    // ------------- IMPLEMENT HERE ---------
    // Compute the normals for each vertex v in the mesh using the constant weights
    // technique (see .pdf) and store it inside v_cste_weights_n[v]
    // ------------- IMPLEMENT HERE ---------

    //allocation of iterating vertices, begin, end (cf slide 35 SM tuto)
    Mesh::Vertex_iterator v_it, v_begin, v_end;

    //init vertices
    v_begin = mesh.vertices_begin();
    v_end = mesh.vertices_end();

    //iterate over mesh vertices
    for( v_it = v_begin ; v_it != v_end; ++ v_it )
    {

        //allocate the circulator
        Mesh::Vertex v = *v_it;
        Mesh::Vertex_around_vertex_circulator vc, vc_end;

        Point pc, p1, p2, crossprod, sum;

        //begin and end
        vc = mesh.vertices(v);
        vc_end = vc;
        //center position
        pc = mesh.position(v);
        sum = (0.0, 0.0, 0.0);

        //loop around vertex v_it (cf slide 40 SM tuto)
        do{
            //p1 and p2 are both neighbours ready for cross product
            p1 = mesh.position(*vc);
            p2 = mesh.position(*(++vc));
            crossprod = cross(p1-pc,p2-pc);
            sum = sum + crossprod/norm(crossprod);
        }while(vc != vc_end);
        v_cste_weights_n[v] = sum/norm(sum);
    }
}
// ========================================================================
// EXERCISE 1.2
// ========================================================================
void Viewer::computeNormalsByAreaWeights() {
    Point default_normal(0.0, 1.0, 0.0);
    Mesh::Vertex_property<Point> v_area_weights_n =
            mesh.vertex_property<Point>("v:area_weight_n", default_normal);

    // ------------- IMPLEMENT HERE ---------
    // Compute the normals for each vertex v in the mesh using the weights proportionals
    // to the areas technique (see .pdf) and store inside v_area_weights_n[v]
    // ------------- IMPLEMENT HERE ---------

    //allocation of iterating vertices, begin, end (cf slide 35 SM tuto)
    Mesh::Vertex_iterator v_it, v_begin, v_end;

    //init vertices
    v_begin = mesh.vertices_begin();
    v_end = mesh.vertices_end();

    //iterate over mesh vertices
    for( v_it = v_begin ; v_it != v_end; ++ v_it )
    {

        //allocation du circulateur
        Mesh::Vertex v = *v_it;
        Mesh::Vertex_around_vertex_circulator vc, vc_end;
        Point pc, p1, p2, crossprod, sum;
        float area = 0;

        //début et fin
        vc = mesh.vertices(v);
        vc_end = vc;
        //position du centre
        pc = mesh.position(v);
        sum = (0.0, 0.0, 0.0);

        //loop around vertex v_it (cf slide 40 SM tuto)
        do{
            //on stock deux voisins à chaque étapes pour le cross product
            p1 = mesh.position(*vc);
            p2 = mesh.position(*(++vc));
            crossprod = cross(p1-pc,p2-pc);
            area = norm(crossprod)/2; //c'est comme ça l'air non?
            sum = sum + crossprod*area/norm(crossprod);
            //version simplifiée sans area :
            //sum = sum + crossprod/2;
        }while(vc != vc_end);
        v_area_weights_n[v] = sum/norm(sum);
    }
}
// ========================================================================
// EXERCISE 1.3
// ========================================================================
void Viewer::computeNormalsWithAngleWeights() {
    Point default_normal(0.0, 1.0, 0.0);
    Mesh::Vertex_property<Point> v_angle_weights_n =
            mesh.vertex_property<Point>("v:angle_weight_n", default_normal);

    // ------------- IMPLEMENT HERE ---------
    // Compute the normals for each vertex v in the mesh using the weights proportionals
    // to the angles technique (see .pdf) and store it inside v_angle_weights_n[v]
    // ------------- IMPLEMENT HERE ---------

    //Iterator vertices allocation
    Mesh::Vertex_iterator v_it, v_begin, v_end;

    //init vertices
    v_begin = mesh.vertices_begin();
    v_end = mesh.vertices_end();

    //Loop over each vertex of the mesh
    for( v_it = v_begin ; v_it != v_end; ++ v_it )
    {

        //allocate the circulator
        Mesh::Vertex v = *v_it;
        Mesh::Vertex_around_vertex_circulator vc, vc_end;
        Point pc, p1, p2, crossprod, sum;
        Scalar angle = 0;
        Scalar dotprod = 0;

        //definition of begin and end points
        vc = mesh.vertices(v);
        vc_end = vc;
        //position of the center
        pc = mesh.position(v);
        sum = (0.0, 0.0, 0.0);

        //circulation loop over the center
        do{
            //Here we define auxiliary variable for more lisibility

            //Saving two neighboors for the crossporduct
            p1 = mesh.position(*vc);
            p2 = mesh.position(*(++vc));

            //crossproduct
            crossprod = cross(p1-pc,p2-pc);

            //dotproduct
            dotprod = dot(p1-pc,p2-pc);
            //from a*b = |a||b|cos(theta)
            angle = acos(dotprod/(norm(p1-pc)*norm(p2-pc)));

            sum = sum + angle * crossprod/norm(crossprod);
        }while(vc != vc_end);
        v_angle_weights_n[v] = sum/norm(sum);
    }
}

// ========================================================================
// EXERCISE 2.1
// ========================================================================

void Viewer::calc_uniform_laplacian() {
    Mesh::Vertex_property<Scalar> v_uniLaplace = mesh.vertex_property<Scalar>("v:uniLaplace", 0);
    Mesh::Vertex_around_vertex_circulator   vv_c, vv_end;
    Point             laplace(0.0);
    min_uniLaplace = 1000;
    max_uniLaplace = -1000;

    // ------------- IMPLEMENT HERE ---------
    // For each non-boundary vertex, compute uniform Laplacian operator vector
    // and store the vector length in the vertex property of the
    // mesh called v_uniLaplace[v].
    // Store min and max values of v_uniLaplace[v] in min_uniLaplace and max_uniLaplace.
    // ------------- IMPLEMENT HERE ---------

    //allocation of iterating vertices, begin, end (cf slide 35 SM tuto)
    Mesh::Vertex_iterator v_it, v_begin, v_end;

    //init vertices
    v_begin = mesh.vertices_begin();
    v_end = mesh.vertices_end();

    //iterate over mesh vertices
    for(v_it = v_begin ; v_it != v_end; ++ v_it )
    {
        //allocation de l'itérateur pour tourner autour de chaque vertex
        Mesh::Vertex v = *v_it;
        //début et fin
        vv_c = mesh.vertices(v);
        vv_end = vv_c;

        Scalar n = 0 ;
        laplace = 0;

        //loop around vertex v_it (cf slide 40 SM tuto)
        do {
            Mesh::Vertex vi = *vv_c;
            laplace += (mesh.position(vi) - mesh.position(v));
            n++;
        } while(++vv_c != vv_end);
        laplace = laplace/n ;

        // Check min or max ?
        if(norm(laplace) > max_uniLaplace){
            max_uniLaplace = norm(laplace) ;
        }
        else if(norm(laplace) < min_uniLaplace){
            min_uniLaplace = norm(laplace) ;
        }
        v_uniLaplace[v] = norm(laplace) ;
    }
}


// ========================================================================

// EXERCISE 2.2

// ========================================================================

void Viewer::calc_mean_curvature() {
    Mesh::Vertex_property<Scalar>  v_curvature = mesh.vertex_property<Scalar>("v:curvature", 0);
    Mesh::Edge_property<Scalar> e_weight = mesh.edge_property<Scalar>("e:weight", 0);
    Mesh::Vertex_property<Scalar>  v_weight = mesh.vertex_property<Scalar>("v:weight", 0);
    Mesh::Halfedge_around_vertex_circulator vh_c, vh_end;
    Mesh::Vertex neighbor_v;
    Mesh::Edge e;
    Point laplace(0.0);
    min_mean_curvature = 1000;
    max_mean_curvature = -1000;

    // ------------- IMPLEMENT HERE ---------
    // For all non-boundary vertices, approximate the mean curvature using
    // the length of the Laplace-Beltrami approximation.
    // Save your approximation in v_curvature vertex property of the mesh.
    // Use the weights from calc_weights(): e_weight and v_weight
    // ------------- IMPLEMENT HERE ---------

    //allocation of iterating vertices, begin, end (cf slide 35 SM tuto)
    Mesh::Vertex_iterator v_it, v_begin, v_end;

    //init vertices
    v_begin = mesh.vertices_begin();
    v_end = mesh.vertices_end();

    //iterate over mesh vertices
    for( v_it = v_begin ; v_it != v_end; ++ v_it )
    {
        //allocate itérator to turn around each vertex
        Mesh::Vertex v = *v_it;
        Point P = mesh.position(v);

        //begin and end
        vh_c = mesh.halfedges(v);
        vh_end = vh_c;
        laplace = 0;

        //loop around vertex v_it (cf slide 40 SM tuto)
        do{
            Point Pi = mesh.position(mesh.to_vertex(*vh_c));
            e = mesh.edge(*vh_c);
            laplace = laplace + (Pi-P)*e_weight[e];
        }while(++vh_c != vh_end);

        laplace = laplace * v_weight[*v_it] ;

        if(norm(laplace) > max_mean_curvature){
            max_mean_curvature = norm(laplace) ;
        }
        else if(norm(laplace) < min_mean_curvature){
            min_mean_curvature = norm(laplace) ;
        }
        v_curvature[*v_it] = norm(laplace) ;
    }
}
// ========================================================================
// EXERCISE 2.3
// ========================================================================

void Viewer::calc_gauss_curvature() {
    Mesh::Vertex_property<Scalar> v_gauss_curvature = mesh.vertex_property<Scalar>("v:gauss_curvature", 0);
    Mesh::Vertex_property<Scalar> v_weight = mesh.vertex_property<Scalar>("v:weight", 0);
    Mesh::Vertex_around_vertex_circulator vv_c, vv_c2, vv_end;
    Point d0, d1;
    Scalar angles, cos_angle;
    Scalar lb(-1.0f), ub(1.0f);
    min_gauss_curvature = 1000;
    max_gauss_curvature = -1000;

    // ------------- IMPLEMENT HERE ---------
    // For each non-boundary vertex, approximate Gaussian curvature,
    // and store it in the vertex property v_gauss_curvature.
    // Hint: When calculating angles out of cross products make sure the value
    // you pass to the acos function is between -1.0 and 1.0.
    // Use the v_weight property for the area weight.
    // ------------- IMPLEMENT HERE ---------

    //allocation of iterating vertices, begin, end (cf slide 35 SM tuto)
    Mesh::Vertex_iterator v_it, v_begin, v_end;

    //init vertices
    v_begin = mesh.vertices_begin();
    v_end = mesh.vertices_end();

    Scalar G = 0;

    //iterate over mesh vertices
    for( v_it = v_begin ; v_it != v_end; ++ v_it )
    {

        //allocate itérator to turn around each vertex
        Point d = mesh.position(*v_it) ;

        //begin and end
        Mesh::Vertex v = *v_it;
        vv_c = mesh.vertices(v);

        vv_end = vv_c;

        angles = 0;

        //loop around vertex v_it (cf slide 40 SM tuto)
        do{
            d0 = mesh.position(*vv_c) ;
            vv_c2 = ++vv_c;
            d1 = mesh.position(*vv_c2) ;
            cos_angle = dot(d-d0,d-d1)/(norm(d-d0) * norm(d-d1)) ;

            //checking if cos_angle is between -1 and 1 and restraining it
            if(cos_angle>ub){ cos_angle = ub; }
            if(cos_angle<lb){ cos_angle = lb; }

            angles += acos(cos_angle);
        }while(vv_c != vv_end);

        G = (2*M_PI -angles)*2*v_weight[*v_it] ;

        if(G > max_gauss_curvature){
            max_gauss_curvature = G ;
        }
        else if(G < min_gauss_curvature){
            min_gauss_curvature = G ;
        }
        v_gauss_curvature[*v_it] = G ;
    }
}
