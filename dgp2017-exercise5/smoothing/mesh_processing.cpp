//=============================================================================
//
//   Code framework for the lecture
//
//   "Digital 3D Geometry Processing"
//
//   Gaspard Zoss, Alexandru Ichim
//
//   Copyright (C) 2016 by Computer Graphics and Geometry Laboratory,
//         EPF Lausanne
//
//   Edited 2017
//-----------------------------------------------------------------------------
#include "mesh_processing.h"
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
// ======================================================================
// EXERCISE 1.1
// ========================================================================
void MeshProcessing::uniform_smooth(const unsigned int iterations) {

    for (unsigned int iter=0; iter<iterations; ++iter) {
        // ------------- IMPLEMENT HERE ---------
        // For each non-boundary vertex, update its position according to the uniform Laplacian operator
        // ------------- IMPLEMENT HERE ---------
        Scalar dtl = 0.2;
        //allocation of iterating vertices, begin, end (cf slide 35 SM tuto)
        Mesh::Vertex_iterator v_it, v_begin, v_end;
        Mesh::Vertex_around_vertex_circulator   vv_c, vv_end;
        Point             laplace(0.0);
        Scalar              n;

        //allocation of a temporary array to stock all p' positions
        Mesh::Vertex_property<Point> points = mesh_.get_vertex_property<Point>("v:point");

        //init vertices
        v_begin = mesh_.vertices_begin();
        v_end = mesh_.vertices_end();
        //iterate over mesh vertices to compute p'
        for(v_it = v_begin ; v_it != v_end; ++ v_it )
        {
            if (!mesh_.is_boundary(*v_it)){
                vv_c = mesh_.vertices(*v_it);
                vv_end = vv_c;
                n = 0;
                do {
                    laplace += (mesh_.position(*vv_c) - mesh_.position(*v_it));
                    ++n;
                } while(++vv_c != vv_end);

                laplace /= n;
                points[*v_it] = mesh_.position(*v_it) + dtl * laplace;
            }
        }
        //iterating a second time to put new p' in p
        for(v_it = v_begin ; v_it != v_end; ++ v_it ){
            mesh_.position(*v_it) = points[*v_it];
        }
    }
}

// ======================================================================
// EXERCISE 1.2
// ========================================================================
void MeshProcessing::smooth(const unsigned int iterations)
{
    for(unsigned int iter=0; iter<iterations; ++iter)
    {
        // ------------- IMPLEMENT HERE ---------
        // Perform Cotan Laplacian smoothing:
        // 1) precompute edge weights using calc_edge_weights()
        // 2) for each non-boundary vertex, update its position using the normalized cotan Laplacian operator
        // (Hint: use the precomputed edge weights in the edge property "e:weight")
        // ------------- IMPLEMENT HERE ---------


        //allocation of iterating vertices, begin, end (cf slide 35 SM tuto) and halfedge circulator
        Mesh::Vertex_iterator v_it, v_begin, v_end;
        Mesh::Halfedge_around_vertex_circulator vv_c, vv_end;

        // Some variable allocation
        Point laplace(0.0);
        Scalar sum_weights;
        Scalar dtl = 0.2;
        Mesh::Edge e;
        Mesh::Vertex ver ;

        //allocation of a temporary array to stock all p' positions
        Mesh::Vertex_property<Point> points = mesh_.get_vertex_property<Point>("v:point");

        //precomputing and stocking the edge weights {w_j}
        calc_edges_weights();
        Mesh::Edge_property<Scalar> weights = mesh_.get_edge_property<Scalar>("e:weight");

        // init vertices
        v_begin = mesh_.vertices_begin();
        v_end = mesh_.vertices_end();

        //iterate over mesh vertices to compute p'
        for(v_it=v_begin ; v_it!=v_end; ++v_it ){
            if (!mesh_.is_boundary(*v_it)){ // boundary check

                // Circulator initialisation
                vv_c = mesh_.halfedges(*v_it);
                vv_end = vv_c;

                // resetting the sums
                laplace = {0,0,0};
                sum_weights = 0;
                do{
                    ver = mesh_.to_vertex(*vv_c) ;
                    e = mesh_.edge(*vv_c);
                    laplace += weights[e]*(mesh_.position(ver) - mesh_.position(*v_it));
                    sum_weights += weights[e];
                } while(++vv_c != vv_end);
                laplace = laplace/sum_weights;
                points[*v_it] = mesh_.position(*v_it) + dtl * laplace;
            }
        }

        //iterating a second time to put new p' in p
        for(v_it = v_begin ; v_it != v_end; ++ v_it ){
            mesh_.position(*v_it) = points[*v_it];
        }
    }
}

// ======================================================================
// EXERCISE 2
// ========================================================================
void MeshProcessing::implicit_smoothing(const double timestep) {

    const int n = mesh_.n_vertices();

    // get vertex position
    auto points = mesh_.vertex_property<Point>("v:point");

    // compute cotan edge weights and vertex areas
    calc_weights ();
    auto cotan = mesh_.edge_property<Scalar>("e:weight");
    auto area_inv = mesh_.vertex_property<Scalar>("v:weight");

    // A*X = B
    Eigen::SparseMatrix<double> A(n,n);
    Eigen::MatrixXd B(n,3);

    // nonzero elements of A as triplets: (row, column, value)
    std::vector< Eigen::Triplet<double> > triplets;

    // ========================================================================
    // TODO: IMPLEMENTATION FOR EXERCISE 2 HERE

    //allocation of iterating vertices, begin, end (cf slide 35 SM tuto)
    Mesh::Vertex_iterator v_it, v_begin, v_end;

    // some variable declarations
    Scalar dtl = 0.2;
    Scalar sum = 0.;
    Mesh::Edge e;
    Mesh::Vertex ver;
    Point p;

    //init vertices
    v_begin = mesh_.vertices_begin();
    v_end = mesh_.vertices_end();

    //iterate over all mesh vertices in order to fill the matrices
    Mesh::Vertex v;
    Mesh::Halfedge_around_vertex_circulator vc, vc_end;
    for(v_it = v_begin; v_it != v_end; ++ v_it) {

        //allocate the halfedge circulator
        v = *v_it;
        sum = 0.;


        //begin and end
        vc = mesh_.halfedges(v);
        vc_end = vc;

        //center position
        p = mesh_.position(v);

        // Fill matrix B (if i understood it correctly)
        B(v.idx(), 0) = p[0] / area_inv[v];
        B(v.idx(), 1) = p[1] / area_inv[v];
        B(v.idx(), 2) = p[2] / area_inv[v];

        //loop around vertex v_it (cf slide 40 SM tuto) to get the necessary values to fill 'M'
        do {
            e = mesh_.edge(*vc);
            ver = mesh_.to_vertex(*vc) ;
            sum += cotan[e] ; // Sum used to obtain M(i,i)
            triplets.push_back(Eigen::Triplet<double>(v.idx(), ver.idx(), (2 / area_inv[v]) - dtl * cotan[e]));
        } while(++vc != vc_end);
        triplets.push_back(Eigen::Triplet<double>(v.idx(), v.idx(), (2 / area_inv[v]) + dtl * sum));
    }

    // ========================================================================

    // build sparse matrix from triplets
    A.setFromTriplets(triplets.begin(), triplets.end());

    // solve A*X = B
    Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > solver(A);
    Eigen::MatrixXd X = solver.solve(B);

    // copy solution
    for (int i = 0; i < n; ++i)
    {
        Mesh::Vertex v(i);
        for (int dim = 0; dim < 3; ++dim)
            points[v][dim] = X(i, dim);
    }

    // clean-up
    mesh_.remove_vertex_property(area_inv);
    mesh_.remove_edge_property(cotan);
}

// ======================================================================
// EXERCISE 3.1
// ========================================================================
void MeshProcessing::uniform_laplacian_enhance_feature(const unsigned int iterations,
                                                       const unsigned int coefficient) {

    // ------------- IMPLEMENT HERE ---------
    // Feature enhancement using the uniform Laplacian operator:
    // 1) perform uniform Laplacian smoothing for enhancement_smoothing_iterations iterations
    // 2) update the vertex positions according to the difference between the original and the smoothed mesh,
    //    using enhancement_coef as the value of alpha in the feature enhancement formula
    // ------------- IMPLEMENT HERE ---------

    Mesh::Vertex_iterator v_it, v_begin, v_end;

    //init vertices
    v_begin = mesh_.vertices_begin();
    v_end = mesh_.vertices_end();

    // get vertex position before smooting
    Point* p_in = new Point[mesh_.n_vertices()];
    int i = 0;
    for (v_it = v_begin; v_it != v_end; ++ v_it) {
        Mesh::Vertex v = *v_it;
        p_in[i] = mesh_.position(v);
        i++;
    }

    // Perform the uniform smoothing
    MeshProcessing::uniform_smooth(iterations);


    // get vertex position after smooting
    Point p_out;
    i = 0;
    for (v_it = v_begin; v_it != v_end; ++ v_it) {
        Mesh::Vertex v = *v_it;
        p_out = mesh_.position(v) ;
        mesh_.position(v) = p_out + coefficient * (p_in[i] - p_out);
        i++;
    }
}

// ======================================================================
// EXERCISE 3.2
// ========================================================================
void MeshProcessing::cotan_laplacian_enhance_feature(const unsigned int iterations,
                                                     const unsigned int coefficient) {

    // ------------- IMPLEMENT HERE ---------
    // Feature enhancement using the normalized cotan Laplacian operator:
    // 1) perform cotan Laplacian smoothing for enhancement_smoothing_iterations iterations
    // 2) update the vertex positions according to the difference between the original and the smoothed mesh,
    //    using enhancement_coef as the value of alpha in the feature enhancement formula
    // ------------- IMPLEMENT HERE ---------

    Mesh::Vertex_iterator v_it, v_begin, v_end;

    //init vertices
    v_begin = mesh_.vertices_begin();
    v_end = mesh_.vertices_end();

    // get vertex position before smooting
    Point* p_in = new Point[mesh_.n_vertices()];
    int i = 0;
    for (v_it = v_begin; v_it != v_end; ++ v_it) {
        Mesh::Vertex v = *v_it;
        p_in[i] = mesh_.position(v);
        i++;
    }

    // Perform the uniform smoothing
    MeshProcessing::smooth(iterations);


    // get vertex position after smooting
    Point p_out;
    i = 0;
    for (v_it = v_begin; v_it != v_end; ++ v_it) {
        Mesh::Vertex v = *v_it;
        p_out = mesh_.position(v) ;
        mesh_.position(v) = p_out + coefficient * (p_in[i] - p_out);
        i++;
    }
}

void MeshProcessing::calc_weights() {
    calc_edges_weights();
    calc_vertices_weights();
}

void MeshProcessing::calc_uniform_mean_curvature() {
    Mesh::Vertex_property<Scalar> v_unicurvature =
            mesh_.vertex_property<Scalar>("v:unicurvature", 0.0f);
    Mesh::Vertex_around_vertex_circulator   vv_c, vv_end;
    Point             laplace(0.0);

    for (auto v: mesh_.vertices()) {
        Scalar curv = 0;

        if (!mesh_.is_boundary(v)) {
            laplace = Point(0.0f);
            double n = 0;
            vv_c = mesh_.vertices(v);
            vv_end = vv_c;

            do {
                laplace += (mesh_.position(*vv_c) - mesh_.position(v));
                ++n;
            } while(++vv_c != vv_end);

            laplace /= n;

            curv = norm(laplace);
        }
        v_unicurvature[v] = curv;
    }
}

void MeshProcessing::calc_mean_curvature() {
    Mesh::Vertex_property<Scalar>  v_curvature =
            mesh_.vertex_property<Scalar>("v:curvature", 0.0f);
    Mesh::Edge_property<Scalar> e_weight =
            mesh_.edge_property<Scalar>("e:weight", 0.0f);
    Mesh::Vertex_property<Scalar>  v_weight =
            mesh_.vertex_property<Scalar>("v:weight", 0.0f);

    Mesh::Halfedge_around_vertex_circulator vh_c, vh_end;
    Mesh::Vertex neighbor_v;
    Mesh::Edge e;
    Point laplace(0.0f, 0.0f, 0.0f);

    for (auto v: mesh_.vertices()) {
        Scalar curv = 0.0f;

        if (!mesh_.is_boundary(v)) {
            laplace = Point(0.0f, 0.0f, 0.0f);

            vh_c = mesh_.halfedges(v);
            vh_end = vh_c;

            do {
                e = mesh_.edge(*vh_c);
                neighbor_v = mesh_.to_vertex(*vh_c);
                laplace += e_weight[e] * (mesh_.position(neighbor_v) -
                                          mesh_.position(v));

            } while(++vh_c != vh_end);

            laplace *= v_weight[v];
            curv = norm(laplace);
        }
        v_curvature[v] = curv;
    }
}

void MeshProcessing::calc_gauss_curvature() {
    Mesh::Vertex_property<Scalar> v_gauss_curvature =
            mesh_.vertex_property<Scalar>("v:gauss_curvature", 0.0f);
    Mesh::Vertex_property<Scalar> v_weight =
            mesh_.vertex_property<Scalar>("v:weight", 0.0f);
    Mesh::Vertex_around_vertex_circulator vv_c, vv_c2, vv_end;
    Point d0, d1;
    Scalar angles, cos_angle;
    Scalar lb(-1.0f), ub(1.0f);

    // compute for all non-boundary vertices
    for (auto v: mesh_.vertices()) {
        Scalar curv = 0.0f;

        if (!mesh_.is_boundary(v)) {
            angles = 0.0f;

            vv_c = mesh_.vertices(v);
            vv_end = vv_c;

            do {
                vv_c2 = vv_c;
                ++ vv_c2;
                d0 = normalize(mesh_.position(*vv_c) - mesh_.position(v));
                d1 = normalize(mesh_.position(*vv_c2) - mesh_.position(v));
                cos_angle = max(lb, min(ub, dot(d0, d1)));
                angles += acos(cos_angle);
            } while(++vv_c != vv_end);

            curv = (2 * 3.1415 - angles) * 2.0f * v_weight[v];
        }
        v_gauss_curvature[v] = curv;
    }
}

void MeshProcessing::calc_edges_weights() {
    auto e_weight = mesh_.edge_property<Scalar>("e:weight", 0.0f);
    auto points = mesh_.vertex_property<Point>("v:point");

    Mesh::Halfedge h0, h1, h2;
    Point p0, p1, p2, d0, d1;

    for (auto e: mesh_.edges())
    {
        e_weight[e] = 0.0;

        h0 = mesh_.halfedge(e, 0);
        p0 = points[mesh_.to_vertex(h0)];

        h1 = mesh_.halfedge(e, 1);
        p1 = points[mesh_.to_vertex(h1)];

        if (!mesh_.is_boundary(h0))
        {
            h2 = mesh_.next_halfedge(h0);
            p2 = points[mesh_.to_vertex(h2)];
            d0 = p0 - p2;
            d1 = p1 - p2;
            e_weight[e] += dot(d0,d1) / norm(cross(d0,d1));
        }

        if (!mesh_.is_boundary(h1))
        {
            h2 = mesh_.next_halfedge(h1);
            p2 = points[mesh_.to_vertex(h2)];
            d0 = p0 - p2;
            d1 = p1 - p2;
            e_weight[e] += dot(d0,d1) / norm(cross(d0,d1));
        }
    }
}

void MeshProcessing::calc_vertices_weights() {
    Mesh::Face_around_vertex_circulator vf_c, vf_end;
    Mesh::Vertex_around_face_circulator fv_c;
    Scalar area;
    auto v_weight = mesh_.vertex_property<Scalar>("v:weight", 0.0f);

    for (auto v: mesh_.vertices()) {
        area = 0.0;
        vf_c = mesh_.faces(v);

        if(!vf_c) {
            continue;
        }

        vf_end = vf_c;

        do {
            fv_c = mesh_.vertices(*vf_c);

            const Point& P = mesh_.position(*fv_c);  ++fv_c;
            const Point& Q = mesh_.position(*fv_c);  ++fv_c;
            const Point& R = mesh_.position(*fv_c);

            area += norm(cross(Q-P, R-P)) * 0.5f * 0.3333f;

        } while(++vf_c != vf_end);

        v_weight[v] = 0.5 / area;
    }
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

    calc_weights();
    calc_uniform_mean_curvature();
    calc_mean_curvature();
    calc_gauss_curvature();
    color_coding(vertex_valence, &mesh_, v_color_valence, 100 /* bound */);
    color_coding(v_unicurvature, &mesh_, v_color_unicurvature);
    color_coding(v_curvature, &mesh_, v_color_curvature);
    color_coding(v_gauss_curvature, &mesh_, v_color_gaussian_curv);

    // get the mesh attributes and upload them to the GPU
    int j = 0;
    unsigned int n_vertices(mesh_.n_vertices());

    // Create big matrices to send the data to the GPU with the required
    // format
    color_valence_ = Eigen::MatrixXf(3, n_vertices);
    color_unicurvature_ = Eigen::MatrixXf(3, n_vertices);
    color_curvature_ = Eigen::MatrixXf(3, n_vertices);
    color_gaussian_curv_ = Eigen::MatrixXf(3, n_vertices);
    normals_ = Eigen::MatrixXf(3, n_vertices);
    points_ = Eigen::MatrixXf(3, n_vertices);
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
        ++j;
    }
}

void MeshProcessing::color_coding(Mesh::Vertex_property<Scalar> prop, Mesh *mesh,
                                  Mesh::Vertex_property<Color> color_prop, int bound) {
    // Get the value array
    std::vector<Scalar> values = prop.vector();

    // discard upper and lower bound
    unsigned int n = values.size()-1;
    unsigned int i = n / bound;
    std::sort(values.begin(), values.end());
    Scalar min_value = values[i], max_value = values[n-1-i];

    // map values to colors
    for (auto v: mesh->vertices())
    {
        set_color(v, value_to_color(prop[v], min_value, max_value), color_prop);
    }
}

void MeshProcessing::set_color(Mesh::Vertex v, const Color& col,
                               Mesh::Vertex_property<Color> color_prop)
{
    color_prop[v] = col;
}

Color MeshProcessing::value_to_color(Scalar value, Scalar min_value, Scalar max_value) {
    Scalar v0, v1, v2, v3, v4;
    v0 = min_value + 0.0/4.0 * (max_value - min_value);
    v1 = min_value + 1.0/4.0 * (max_value - min_value);
    v2 = min_value + 2.0/4.0 * (max_value - min_value);
    v3 = min_value + 3.0/4.0 * (max_value - min_value);
    v4 = min_value + 4.0/4.0 * (max_value - min_value);

    Color col(1.0f, 1.0f, 1.0f);

    if (value < v0) {
        col = Color(0, 0, 1);
    } else if (value > v4) {
        col = Color(1, 0, 0);
    } else if (value <= v2) {
        if (value <= v1) { // [v0, v1]
            Scalar u =  (value - v0) / (v1 - v0);
            col = Color(0, u, 1);
        } else { // ]v1, v2]
            Scalar u = (value - v1) / (v2 - v1);
            col = Color(0, 1, 1-u);
        }
    } else {
        if (value <= v3) { // ]v2, v3]
            Scalar u = (value - v2) / (v3 - v2);
            col = Color(u, 1, 0);
        } else { // ]v3, v4]
            Scalar u = (value - v3) / (v4 - v3);
            col = Color(1, 1-u, 0);
        }
    }
    return col;
}

MeshProcessing::~MeshProcessing() {}
}
