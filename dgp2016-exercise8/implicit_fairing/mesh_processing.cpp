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
#define _USE_MATH_DEFINES
#include "mesh_processing.h"
#include <cmath>
#include <set>

namespace mesh_processing {

using surface_mesh::Point;
using surface_mesh::Scalar;
using surface_mesh::Color;
using std::min;
using std::max;
using std::cout;
using std::endl;

MeshProcessing::MeshProcessing(const string& filename) {
    load_mesh(filename);
}

void MeshProcessing::create_dual_mesh() { //"DONE"

    Mesh original_mesh_ = mesh_;
    mesh_.clear();

    Mesh::Face_around_vertex_circulator  vf_c;
    Mesh::Face_around_vertex_circulator vf_end;
    Mesh::Vertex_around_face_circulator fv_c;
    Scalar it = 0 ;

    for (Mesh::Vertex_iterator v_it = original_mesh_.vertices_begin(); v_it != original_mesh_.vertices_end(); ++v_it) {

        // ------------- IMPLEMENT HERE ---------
        // For each vertex:
        //     compute centers of mass of the faces around vertex
        //     add them to the mesh as vertices
        //     create a face that consists of all these vertices that are centers of the faces around the current vertex
        // ------------- IMPLEMENT HERE ---------
        Mesh::Vertex_around_face_circulator fv_end;
        //The future new vector in the mesh.
        Mesh::Vertex new_v;
        //The vector containing the vertices for the new face
        std::vector<Mesh::Vertex> new_face_vertices;
        //center of mass
        Point cm;
        //counter of vertices arround a face
        Scalar Nv = 0;

        //Inititalizing the loop over faces arround vertex v_it
        vf_c = original_mesh_.faces(*v_it);
        vf_end = vf_c;
        //Starting the loop over faces arround vertex v_it
        do{
            Nv = 0;
            cm = 0;

            //Inititalizing the loop over vertices arround face vf_c
            fv_c = original_mesh_.vertices(*vf_c);
            fv_end = fv_c;
            //Starting the loop over vertices arround face vf_c
            do{
                cm += original_mesh_.position(*fv_c);
                Nv++;
            }while(++fv_c != fv_end);
            ++it;
            cout << it << endl ;
            //finishing the average of vertices positions
            cm /= Nv;
            //adding the new position as a new vertex on mesh_ (which is cleared in the beginning)
            new_v = mesh_.add_vertex(cm);
            //Adding the new vertex to the face vertices vector.
            new_face_vertices.push_back(new_v);

        } while(++vf_c != vf_end);

        //adding the new face to the mesh.
        mesh_.add_face(new_face_vertices);
        new_face_vertices.clear();
    }

}

void MeshProcessing::add_primitive(const Eigen::Matrix3f & R, const Eigen::Vector3f & t, const Mesh & mesh_primitive_) {//"DONE"
    Mesh::Vertex_around_face_circulator fv_c;
    Mesh::Vertex_around_face_circulator fv_end;

    for (Mesh::Face_iterator f_it = mesh_primitive_.faces_begin(); f_it != mesh_primitive_.faces_end(); ++f_it) {

        // ------------- IMPLEMENT HERE ---------
        // For each face of cylinder/sphere mesh:
        //	   For each vertex of the current face:
        //          apply the rotation R and translation t to the vertex position
        //          add the vertex to the mesh
        //     add the face that consists of the new vertices to the mesh
        // ------------- IMPLEMENT HERE ---------

        Eigen::Vector3f Eigenp;
        Point p;
        std::vector<Mesh::Vertex> new_face;
        Mesh::Vertex v;

        //Inititalizing the loop over vertices arround face vf_c
        fv_c = mesh_primitive_.vertices(*f_it);
        fv_end = fv_c;
        //Starting the loop over vertices arround face vf_c
        do{
            p = mesh_primitive_.position(*fv_c);

            //copying to Eigen vector for computations
            Eigenp(0) = p[0];
            Eigenp(1) = p[1];           // There must be something more efficient to transpose Point to Eigen...
            Eigenp(2) = p[2];

            //applying rotation and translation
            Eigenp = R*Eigenp + t;

            //copying back to Point
            p[0] = Eigenp(0);
            p[1] = Eigenp(1);
            p[2] = Eigenp(2);

            //Adding the position to the mesh
            v = mesh_.add_vertex(p); //Does the mesh automatically delete the old one?...
            //Adding new face's vertices
            new_face.push_back(v);

        }while(++fv_c != fv_end);

        //adding the new face to the mesh.
        mesh_.add_face(new_face);
        new_face.clear();
    }
}

void MeshProcessing::create_cylinder_edges() {

    Mesh mesh_original_ = mesh_;
    mesh_.clear();
    const float CYLINDER_LENGTH = 2;

    float cylinder_radius = 0.015; // adjust for different meshes
    float sphere_radius = 0.02;	// adjust for different meshes

    size_t i = 0;
    size_t num_egdes = mesh_original_.n_edges();
    size_t num_vertices = mesh_original_.n_vertices();

    //ADDED
    Eigen::Matrix3f S;
    Eigen::Matrix3f R;
    Eigen::Vector3f t;
    //Eigen::Vector3f v_d_cylindre;
    Point v_d_cylindre;

    Scalar l;
    Scalar theta;
    Point tp;
    Mesh::Vertex v0;
    Mesh::Vertex v1;
    Point p0;
    Point p1;

    Scalar cos_theta, sin_theta;
    Scalar one_co_th;
    Scalar v0_sqr, v1_sqr, v2_sqr;
    //

    for (Mesh::Edge_iterator e_it = mesh_original_.edges_begin(); e_it != mesh_original_.edges_end(); ++e_it) {

        // ------------- IMPLEMENT HERE ---------
        // For each edge of the original mesh:
        //     Compute scaling matrix S to reduce and stretch the cylinder to the size of an edge
        //	   Compute rotation that aligns the axis of the template cylinder (0, 1, 0) with the edge direction
        //     Compute translation vector for the template cylinder
        //     Call add_primitive function to add the cylinder to the mesh "add_primitive(R * S, t, mesh_cylinder_)"
        // ------------- IMPLEMENT HERE ---------

        cout << "edge " << ++i << " of " << num_egdes << endl;
        // Finding the position of the two vertex of this edge
        v0 = mesh_original_.vertex(*e_it,0);
        v1 = mesh_original_.vertex(*e_it,1);

        p0 = mesh_original_.position(v0) ;
        p1 = mesh_original_.position(v1) ;

        //Computing the scaling matrix S :
        l = mesh_original_.edge_length(*e_it) / CYLINDER_LENGTH;

        S << cylinder_radius, 0, 0,
                0, l, 0,
                0, 0, cylinder_radius;

        // creating direction vector of the cylindre and the edge
        auto v_d = (p0 - p1) / norm(p0 - p1) ;

        v_d_cylindre[0] = 0 ;
        v_d_cylindre[1] = 1 ;
        v_d_cylindre[2] = 0 ;

        //calculating the normal vector
        auto v_perp = cross(v_d, v_d_cylindre) ;
        v_perp = normalize(v_perp);

        // Calculating the angle between the two direction vector
        theta = acos(dot(v_d, v_d_cylindre));

        // Correction angle
        if(theta > M_PI/2){
            theta =  M_PI - theta ;
        }

        // Some variables to avoid code duplication
        cos_theta = cos(theta);
        sin_theta = sin(theta);
        one_co_th = 1 - cos_theta;
        v0_sqr = pow(v_perp[0],2);
        v1_sqr = pow(v_perp[1],2);
        v2_sqr = pow(v_perp[2],2);

        // Rotation matrix around any axe passing through the origin
        R <<    v0_sqr + (1 - v0_sqr) * cos_theta,
                v_perp[0] * v_perp[1] * (one_co_th) - v_perp[2] * sin_theta,
                v_perp[0] * v_perp[2] * (one_co_th) + v_perp[1] * sin_theta,

                v_perp[1] * v_perp[0] * (one_co_th) + v_perp[2] * sin_theta,
                v1_sqr + (1 - v1_sqr) * cos_theta,
                v_perp[1] * v_perp[2] * (one_co_th) - v_perp[0] * sin_theta,

                v_perp[2] * v_perp[0] * (one_co_th) - v_perp[1] * sin_theta,
                v_perp[2] * v_perp[1] * (one_co_th) + v_perp[0] * sin_theta,
                v2_sqr + (1 - v2_sqr) * cos_theta;

        //Computing the translation vector t :
        t(0) = p0[0];
        t(1) = p0[1];
        t(2) = p0[2];

        add_primitive(R * S, t, mesh_cylinder_);
    }

    i = 0;
    for (Mesh::Vertex_iterator v_it = mesh_original_.vertices_begin(); v_it != mesh_original_.vertices_end(); ++v_it) {

        // ------------- IMPLEMENT HERE ---------
        // For each vertex of the original mesh:
        //     Compute scaling matrix S
        //     Compute translation vector for the template sphere
        //     Call add_primitive function to add the sphere to the mesh "add_primitive(S, t, mesh_sphere_)"
        // ------------- IMPLEMENT HERE ---------

        //Computing the scaling matrix S :
        S << sphere_radius, 0, 0,
                0, sphere_radius, 0,
                0, 0, sphere_radius;

        //Computing the translation vector t :
        tp = mesh_original_.position(*v_it);
        t(0) = tp[0];
        t(1) = tp[1];
        t(2) = tp[2];
        add_primitive(S, t, mesh_sphere_);

    }
}

void MeshProcessing::write_mesh() {
    mesh_.write("../data/result.obj");
}

void MeshProcessing::load_mesh(const string &filename) {
    if (!mesh_.read(filename)) {
        std::cerr << "Mesh not found, exiting." << std::endl;
        exit(-1);
    }

    cout << "Mesh "<< filename << " loaded." << endl;
    cout << "# of vertices : " << mesh_.n_vertices() << endl;
    cout << "# of faces : " << mesh_.n_faces() << endl;
    cout << "# of edges : " << mesh_.n_edges() << endl;

    // Compute the center of the mesh
    mesh_center_ = Point(0.0f, 0.0f, 0.0f);
    for (auto v: mesh_.vertices()) {
        mesh_center_ += mesh_.position(v);
    }
    mesh_center_ /= mesh_.n_vertices();

    // Compute the maximum distance from all points in the mesh and the center
    dist_max_ = 0.0f;
    for (auto v: mesh_.vertices()) {
        if (distance(mesh_center_, mesh_.position(v)) > dist_max_) {
            dist_max_ = distance(mesh_center_, mesh_.position(v));
        }
    }

    compute_mesh_properties();

    // Store the original mesh, this might be useful for some computations
    mesh_init_ = mesh_;
}

void MeshProcessing::load_primitives(const string &cylinder_mesh_filename, const string& sphere_mesh_filename) {
    if (!mesh_cylinder_.read(cylinder_mesh_filename)) {
        std::cerr << "Cylinder mesh not found, exiting." << std::endl;
        exit(-1);
    }
    if (!mesh_sphere_.read(sphere_mesh_filename)) {
        std::cerr << "Sphere mesh not found, exiting." << std::endl;
        exit(-1);
    }
}

void MeshProcessing::compute_mesh_properties() {
    Mesh::Vertex_property<Point> vertex_normal =
            mesh_.vertex_property<Point>("v:normal");
    mesh_.update_face_normals();
    mesh_.update_vertex_normals();
    Mesh::Vertex_property<Color> v_color_valence =
            mesh_.vertex_property<Color>("v:color_valence",
                                         Color(1.0f, 1.0f, 1.0f));
    Mesh::Vertex_property<Color> v_color_unicurvature =
            mesh_.vertex_property<Color>("v:color_unicurvature",
                                         Color(1.0f, 1.0f, 1.0f));
    Mesh::Vertex_property<Color> v_color_curvature =
            mesh_.vertex_property<Color>("v:color_curvature",
                                         Color(1.0f, 1.0f, 1.0f));
    Mesh::Vertex_property<Color> v_color_gaussian_curv =
            mesh_.vertex_property<Color>("v:color_gaussian_curv",
                                         Color(1.0f, 1.0f, 1.0f));
    Mesh::Vertex_property<Color> v_color_laplacian =
            mesh_.vertex_property<Color>("v:color_laplacian",
                                         Color(1.0f, 1.0f, 1.0f));

    Mesh::Vertex_property<Scalar> vertex_valence =
            mesh_.vertex_property<Scalar>("v:valence", 0.0f);
    for (auto v: mesh_.vertices()) {
        vertex_valence[v] = mesh_.valence(v);
    }

    Mesh::Vertex_property<Scalar> v_unicurvature =
            mesh_.vertex_property<Scalar>("v:unicurvature", 0.0f);
    Mesh::Vertex_property<Scalar> v_curvature =
            mesh_.vertex_property<Scalar>("v:curvature", 0.0f);
    Mesh::Vertex_property<Scalar> v_gauss_curvature =
            mesh_.vertex_property<Scalar>("v:gauss_curvature", 0.0f);
    Mesh::Vertex_property<Scalar> v_harmonic_function =
            mesh_.vertex_property<Scalar>("v:harmonic_function", 0.0f);

    // get the mesh attributes and upload them to the GPU
    int j = 0;
    unsigned int n_vertices(mesh_.n_vertices());

    // Create big matrices to send the data to the GPU with the required format
    color_valence_ = Eigen::MatrixXf(3, n_vertices);
    color_unicurvature_ = Eigen::MatrixXf(3, n_vertices);
    color_curvature_ = Eigen::MatrixXf(3, n_vertices);
    color_gaussian_curv_ = Eigen::MatrixXf(3, n_vertices);
    color_laplacian_ = Eigen::MatrixXf(3, n_vertices);
    normals_ = Eigen::MatrixXf(3, n_vertices);
    points_ = Eigen::MatrixXf(3, n_vertices);
    selection_ = Eigen::MatrixXf(3, 2);
    indices_ = MatrixXu(3, mesh_.n_faces());

    for(auto f: mesh_.faces()) {
        std::vector<float> vv(3);
        int k = 0;
        for (auto v: mesh_.vertices(f)) {
            vv[k] = v.idx();
            ++k;
        }
        indices_.col(j) << vv[0], vv[1], vv[2];
        ++j;
    }

    j = 0;
    for (auto v: mesh_.vertices()) {
        points_.col(j) << mesh_.position(v).x,
                mesh_.position(v).y,
                mesh_.position(v).z;

        normals_.col(j) << vertex_normal[v].x,
                vertex_normal[v].y,
                vertex_normal[v].z;

        color_valence_.col(j) << v_color_valence[v].x,
                v_color_valence[v].y,
                v_color_valence[v].z;

        color_unicurvature_.col(j) << v_color_unicurvature[v].x,
                v_color_unicurvature[v].y,
                v_color_unicurvature[v].z;

        color_curvature_.col(j) << v_color_curvature[v].x,
                v_color_curvature[v].y,
                v_color_curvature[v].z;

        color_gaussian_curv_.col(j) << v_color_gaussian_curv[v].x,
                v_color_gaussian_curv[v].y,
                v_color_gaussian_curv[v].z;

        color_laplacian_.col(j) << v_color_laplacian[v].x,
                v_color_laplacian[v].y,
                v_color_laplacian[v].z;
        ++j;
    }
}

void MeshProcessing::set_color(Mesh::Vertex v, const Color& col,
                               Mesh::Vertex_property<Color> color_prop)
{
    color_prop[v] = col;
}

MeshProcessing::~MeshProcessing() {}
}
