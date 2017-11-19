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
//-----------------------------------------------------------------------------
#ifndef MESH_PROCESSING_H
#define MESH_PROCESSING_H

#include <surface_mesh/Surface_mesh.h>
#include <Eigen/Sparse>
#include <Eigen/Geometry>

typedef surface_mesh::Surface_mesh Mesh;
typedef Eigen::Matrix<uint32_t, Eigen::Dynamic, Eigen::Dynamic> MatrixXu;

namespace mesh_processing {

using std::string;

class MeshProcessing {

public:
    MeshProcessing(const string& filename);
    ~MeshProcessing();

    const surface_mesh::Point get_mesh_center() { return mesh_center_; }
    const float get_dist_max() { return dist_max_; }
    const Eigen::MatrixXf* get_points() { return &points_; }
	const Eigen::MatrixXf* get_selection() { return &selection_; }
	void set_selection(const Eigen::Vector3f & point, size_t edit_constraint_index) { selection_.col(edit_constraint_index) = point; }
    const MatrixXu* get_indices() { return &indices_; }
    const Eigen::MatrixXf* get_normals() { return &normals_; }
    const Eigen::MatrixXf* get_colors_valence() { return &color_valence_; }
    const Eigen::MatrixXf* get_colors_unicurvature() { return &color_unicurvature_; }
    const Eigen::MatrixXf* get_colors_gaussian_curv() { return &color_gaussian_curv_; }
    const Eigen::MatrixXf* get_color_curvature() { return &color_curvature_; }
	const Eigen::MatrixXf* get_colors_laplacian() { return &color_laplacian_; }
    const unsigned int get_number_of_face() { return mesh_.n_faces(); }
	const unsigned int get_number_of_vertices() { return mesh_.n_vertices(); }

    void load_mesh(const string& filename);
	void write_mesh();
	void load_primitives(const string& cylinder_mesh_filename, const string& sphere_mesh_filename);
    void compute_mesh_properties();

	void create_cylinder_edges();
	void create_dual_mesh();
private:

private:
    Mesh mesh_;
	Mesh mesh_cylinder_;
	Mesh mesh_sphere_;
    Mesh mesh_init_;
    surface_mesh::Point mesh_center_ = surface_mesh::Point(0.0f, 0.0f, 0.0f);
    float dist_max_ = 0.0f;

    Eigen::MatrixXf points_;
	Eigen::MatrixXf selection_;
    MatrixXu indices_;
    Eigen::MatrixXf normals_;
    Eigen::MatrixXf color_valence_;
    Eigen::MatrixXf color_unicurvature_;
    Eigen::MatrixXf color_gaussian_curv_;
    Eigen::MatrixXf color_curvature_;
	Eigen::MatrixXf color_laplacian_;

    void set_color(Mesh::Vertex v, const surface_mesh::Color& col,
                   Mesh::Vertex_property<surface_mesh::Color> color_prop);
	void add_primitive(const Eigen::Matrix3f & T, const Eigen::Vector3f & t, const Mesh & mesh_primitive_);
    };


}

#endif // MESH_PROCESSING_H
