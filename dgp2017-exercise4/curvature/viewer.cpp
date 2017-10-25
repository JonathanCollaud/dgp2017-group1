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

    //allocations des vertices itérateur, début, fin (cf slide 35 SM tuto)
    Mesh::Vertex_iterator v_it, v_begin, v_end;

    //initialisation des vertex alloués ci-dessus
    v_begin = mesh.vertices_begin();
    v_end = mesh.vertices_end();

    //boucle sur chaque vertex du mesh Mesh
    for( v_it = v_begin ; v_it != v_end; ++ v_it )
    {

       //allocation du circulateur
       Mesh::Vertex v = *v_it;
       Mesh::Vertex_around_vertex_circulator vc, vc_end;
       Point pc, p1, p2, crossprod, sum;

       //début et fin
       vc = mesh.vertices(v);
       vc_end = vc;
       //position du centre
       pc = mesh.position(v);
       sum = (0.0, 0.0, 0.0);

       //boucle qui tourne autour du vertex v_it (cf slide 40 SM tuto)
       do{
           //on stock deux voisins à chaque étapes pour le cross product
           p1 = mesh.position(*vc);
           p2 = mesh.position(*(++vc));
           crossprod = cross(p1-pc,p2-pc);
           sum = sum + crossprod/norm(crossprod);
       }while(++vc != vc_end);
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

    //allocations des vertices itérateur, début, fin (cf slide 35 SM tuto)
    Mesh::Vertex_iterator v_it, v_begin, v_end;

    //initialisation des vertex alloués ci-dessus
    v_begin = mesh.vertices_begin();
    v_end = mesh.vertices_end();

    //boucle sur chaque vertex du mesh Mesh
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

       //boucle qui tourne autour du vertex v_it (cf slide 40 SM tuto)
       do{
           //on stock deux voisins à chaque étapes pour le cross product
           p1 = mesh.position(*vc);
           p2 = mesh.position(*(++vc));
           crossprod = cross(p1-pc,p2-pc);
           area = norm(crossprod)/2; //c'est comme ça l'air non?
           sum = sum + crossprod*area/norm(crossprod);
           //version simplifiée sans area :
           //sum = sum + crossprod/2;
       }while(++vc != vc_end);
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

    //allocations des vertices itérateur, début, fin (cf slide 35 SM tuto)
    Mesh::Vertex_iterator v_it, v_begin, v_end;

    //initialisation des vertex alloués ci-dessus
    v_begin = mesh.vertices_begin();
    v_end = mesh.vertices_end();

    //---------------TESTS---------------
    double epsilon = 0.0001; //utilisé pour spot les divisions par 0
    //-----------------------------------

    //boucle sur chaque vertex du mesh Mesh
    for( v_it = v_begin ; v_it != v_end; ++ v_it )
    {

       //allocation du circulateur
       Mesh::Vertex v = *v_it;
       Mesh::Vertex_around_vertex_circulator vc, vc_end;
       Point pc, p1, p2, crossprod, sum;
       float angle = 0;
       float dotprod = 0;

       //début et fin
       vc = mesh.vertices(v);
       vc_end = vc;
       //position du centre
       pc = mesh.position(v);
       sum = (0.0, 0.0, 0.0);

       //boucle qui tourne autour du vertex v_it (cf slide 40 SM tuto)
       do{
           //on stock deux voisins à chaque étapes pour le cross product
           p1 = mesh.position(*vc);
           p2 = mesh.position(*(++vc));
           crossprod = cross(p1-pc,p2-pc);
           dotprod = dot(p1-pc,p2-pc);
           //on ressort :|a^b| = |a||b|sin(theta) -> theta = arcsin(|a^b|/|a||b|)
           //angle = asin(norm(crossprod)/(norm(p1-pc)*norm(p2-pc)));
           //essai avec a*b = |a||b|cos(theta)
           angle = acos(dotprod/(norm(p1-pc)*norm(p2-pc)));

           //-----------------TESTS-------------------
           //on vérifie s'il n'y a pas de division par un nombre trop petit.
           if(norm(p1-pc)<epsilon)
           {cout<< "waring, division by < "<<epsilon<< "because of p1-pc"<<endl;}
           if(norm(p2-pc)<epsilon)
           {cout<< "waring, division by < "<<epsilon<< "because of p2-pc"<<endl;}
           //-----------------------------------------

           sum = sum + crossprod*angle/norm(crossprod);
       }while(++vc != vc_end);
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

    //allocations des vertices itérateur, début, fin (cf slide 35 SM tuto)
        Mesh::Vertex_iterator v_it, v_begin, v_end;

        //initialisation des vertex alloués ci-dessus
        v_begin = mesh.vertices_begin();
        v_end = mesh.vertices_end();

        //boucle sur chaque vertex du mesh Mesh
        for( v_it = v_begin ; v_it != v_end; ++ v_it )
        {
                                                                            // Comment verifier si le vertex est sur une boudary ?


           //allocation de l'itérateur pour tourner autour de chaque vertex
           Mesh::Vertex v = *v_it;
           //Mesh::Vertex_around_vertex_circulator vc, vc_end;

           //début et fin
           vv_c = mesh.vertices(*v_it);
           vv_end = vv_c;

            auto V = mesh.position(v) ;
            Scalar sum = 0;                                              // coment réinitialiser à zéro ?  le type 'vertex' peut ne pas être le bon type ... dépends de la somme
            Scalar n = 0 ;
           //boucle qui tourne autour du vertex v_it (cf slide 40 SM tuto)
           do{
               auto Vi = mesh.position(*vv_c) ;
               n++ ;
               sum = sum + (Vi-V) ;                                         // On est sencé obtenir une valeur scalaire, comment est-ce possible; formule 'v_i - v' qui sont des point en 3D ?!
           }while(++vv_c != vv_end);
            auto Lu = sum/n ;

            // Check min or max ?
            if(Lu > max_uniLaplace){
                max_uniLaplace = Lu ;
            }
            else if(Lu < min_uniLaplace){
                min_uniLaplace = Lu ;
            }
            v_uniLaplace[*v_it] = Lu ;
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
    //Mesh::Vertex neighbor_v;
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


    //allocations des vertices itérateur, début, fin (cf slide 35 SM tuto)
        Mesh::Vertex_iterator v_it, v_begin, v_end;


        //initialisation des vertex alloués ci-dessus
        v_begin = mesh.vertices_begin();
        v_end = mesh.vertices_end();

        //boucle sur chaque vertex du mesh Mesh
        for( v_it = v_begin ; v_it != v_end; ++ v_it )
        {
                                                                            // Comment verifier si le vertex est sur une boudary ?


           //allocation de l'itérateur pour tourner autour de chaque vertex
           auto V = mesh.position(*v_it) ;

           //début et fin
           vh_c = mesh.halfedges(*v_it);
           vh_end = vh_c;

            Scalar sum = 0;
           //boucle qui tourne autour du vertex v_it (cf slide 40 SM tuto)
           do{
               auto vhc = mesh.to_vertex(*vh_c) ;
               auto Vi = mesh.position(vhc) ;
               e = mesh.edge(*vh_c) ;
               sum = sum + (V-Vi)*e_weight[e] ;                             // On est sencé obtenir une valeur scalaire, comment est-ce possible; formule 'v_i - v' qui sont des point en 3D ?!
           }while(++vh_c != vh_end);
            auto Lb = sum*v_weight[*v_it] ;


            if(Lb > max_mean_curvature){
                max_mean_curvature = Lb ;
            }
            else if(Lb < min_mean_curvature){
                min_mean_curvature = Lb ;
            }
            v_curvature[*v_it] = Lb ;
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

    //allocations des vertices itérateur, début, fin (cf slide 35 SM tuto)
        Mesh::Vertex_iterator v_it, v_begin, v_end;


        //initialisation des vertex alloués ci-dessus
        v_begin = mesh.vertices_begin();
        v_end = mesh.vertices_end();

        //boucle sur chaque vertex du mesh Mesh
        for( v_it = v_begin ; v_it != v_end; ++ v_it )
        {
                                                                            // Comment verifier si le vertex est sur une boudary ?


           //allocation de l'itérateur pour tourner autour de chaque vertex

           auto V = mesh.position(*v_it) ;
           //début et fin
           Mesh::Vertex v = *v_it;
           vv_c = mesh.vertices(v);

           vv_end = vv_c;
            Scalar theta_tot = 0;
           //boucle qui tourne autour du vertex v_it (cf slide 40 SM tuto)
           do{
                auto V_1 = mesh.position(*vv_c) ;
                auto V_2 = mesh.position(*++vv_c) ;
                auto vect_1 = V_1 - V ;
                auto vect_2 = V_2 - V ;
                auto theta = asin(norm(cross(vect_1,vect_2))/(norm(vect_1) * norm(vect_2))) ;
                theta_tot = theta_tot + theta ;
                auto G = (2*3.1415 -theta_tot)/v_weight[*v_it] ;


                if(G > max_gauss_curvature){
                    max_gauss_curvature = G ;
                }
                else if(G < min_gauss_curvature){
                    min_gauss_curvature = G ;
                }
            v_gauss_curvature[*v_it] = G ;
            }
        }
}
