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
	Eigen::Vector3f get_closest_vertex(const Eigen::Vector3f & origin, const Eigen::Vector3f & direction, size_t & i);
	int get_closest_face(const Eigen::Vector3f & origin, const Eigen::Vector3f & direction, int mode);
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
    void compute_mesh_properties();

	void deformation();
	void deformation_axis(int mode);
	void calc_mean_curvature();
    void calc_uniform_mean_curvature();
    void calc_gauss_curvature();

	std::vector<surface_mesh::Point> displacement_points_;
	std::vector<surface_mesh::Point> fixed_faces_points_;
	std::vector<size_t> fixed_faces_points_indices_;
	std::vector<surface_mesh::Point> shifted_faces_points_;
	std::vector<size_t> shifted_faces_points_indices_;
	surface_mesh::Point displacement_;

private:
    void calc_weights();
    void calc_edges_weights();
    void calc_vertices_weights();


private:
    Mesh mesh_;
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
	

    void color_coding(Mesh::Vertex_property<surface_mesh::Scalar> prop,
                      Mesh *mesh,
                      Mesh::Vertex_property<surface_mesh::Color> color_prop,
                      int bound = 20);
    void set_color(Mesh::Vertex v, const surface_mesh::Color& col,
                   Mesh::Vertex_property<surface_mesh::Color> color_prop);
    surface_mesh::Color value_to_color(surface_mesh::Scalar value,
                                       surface_mesh::Scalar min_value,
                                       surface_mesh::Scalar max_value);
    };

}

#endif // MESH_PROCESSING_H
