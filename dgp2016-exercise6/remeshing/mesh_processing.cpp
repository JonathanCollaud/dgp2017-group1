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
#include <set>
#include <cmath>
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

	MeshProcessing::~MeshProcessing() {
		// TODO
	}

	void MeshProcessing::remesh(const REMESHING_TYPE &remeshing_type,
		const int &num_iterations) {
		calc_weights();
		calc_mean_curvature();
		calc_uniform_mean_curvature();
		calc_gauss_curvature();
		calc_target_length(remeshing_type);

		// main remeshing loop
		for (int i = 0; i < num_iterations; ++i)
		{
            split_long_edges(); // Works alone, with average and with curvature (do nothing with curvature)
            collapse_short_edges(); // Works alone for average but some bugs for curvature
            equalize_valences(); // The worst one, sry Planck :(
            tangential_relaxation(); // Will it save the rest?
		}
	}

    void MeshProcessing::calc_target_length(const REMESHING_TYPE &remeshing_type) {
        Mesh::Vertex_iterator        v_it,v_begin, v_end(mesh_.vertices_end());
        Mesh::Vertex_around_vertex_circulator  vv_c, vv_end;
        Scalar                   length;
        Scalar                   mean_length;
        Scalar                   H;
        Scalar                   K;

        Mesh::Vertex_property<Scalar> curvature = mesh_.vertex_property<Scalar>("v:meancurvature", 0);
        Mesh::Vertex_property<Scalar> gauss_curvature = mesh_.vertex_property<Scalar>("v:gausscurvature", 0);
        Mesh::Vertex_property<Scalar> target_length = mesh_.vertex_property<Scalar>("v:length", 0);
        Mesh::Vertex_property<Scalar> target_new_length = mesh_.vertex_property<Scalar>("v:newlength", 0);

        const float TARGET_LENGTH = 2.0;

        if (remeshing_type == AVERAGE)
        {
            for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it)
                target_length[*v_it] = TARGET_LENGTH;

        }
        else if (remeshing_type == CURV)
        {
            // ------------- IMPLEMENT HERE ---------
            // Get the maximal curvature at each vertex (use the precomputed mean and gaussian curvature)
            // Calculate the desired edge length as the TARGET_LENGTH divided by the maximal curvature at each vertex, and assign it to the property target_length
            // Smooth the maximal curvature uniformly, use the property vnewtargetlength_ to store the smoothed values intermediately
            // Rescale the property target_new_length such that it's mean equals the user specified TARGET_LENGTH
            // ------------- IMPLEMENT HERE ---------

            Scalar k_max;

            // Calculate mean and gauss curvatures
            MeshProcessing::calc_mean_curvature() ;
            MeshProcessing::calc_gauss_curvature() ;

            //-1-// calculate desired adaptative length
            for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it) {

                if (!mesh_.is_boundary(*v_it)) {
                    // Get the curvatures
                    H = curvature[*v_it] ;
                    K = gauss_curvature[*v_it] ;
                    k_max = H + sqrt(abs(pow(H,2) - K)) ; // Calculate k_max
                    length = TARGET_LENGTH/k_max ; // adapt the target length
                }else
                {
                    length = TARGET_LENGTH; // In case we are on a boundary set the basic target length
                }
                target_length[*v_it] = length;

            }

            //-2-// smooth desired length (used the same scheme as the uniform smoothing performed on mesh but here is performed on lengths)
            Scalar dtl = 0.2;
            Scalar laplace(0.0);
            Scalar n;

            //initialize vertices
            v_begin = mesh_.vertices_begin();
            v_end = mesh_.vertices_end();

            for (int i = 0; i < 5; i++) {
                //iterate over mesh vertices to compute new target_length with the uniform smoothing formula
                for(v_it = v_begin ; v_it != v_end; ++ v_it )
                {
                    if (!mesh_.is_boundary(*v_it)){
                        vv_c = mesh_.vertices(*v_it);
                        vv_end = vv_c;
                        n = 0;
                        laplace = 0 ;
                        do { // sum the difference of all neighbours
                            laplace += (target_length[*vv_c] - target_length[*v_it]);
                            ++n;
                        } while(++vv_c != vv_end);
                        laplace /= n; // get the mean of the difference
                        target_new_length[*v_it] = target_length[*v_it] + dtl * laplace; // Uniform smoothing algorithm
                    }else{
                        target_new_length[*v_it] = target_length[*v_it] ; // on boundary keep the basic target length
                    }
                }
                //Iterating a second time to put new length in target_length
                for(v_it = v_begin ; v_it != v_end; ++ v_it ){
                    target_length[*v_it] = target_new_length[*v_it];
                }
            }

            //-3-// rescale desired length
            Scalar sum = 0 ;
            Scalar ratio ;
            const int N = mesh_.n_vertices();

            for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it) { // loop over all to obtain the mean finally
                sum = sum + target_length[*v_it] ;
            }
            mean_length = sum/N ; // calculate mean length
            ratio = TARGET_LENGTH/mean_length ; // calculate the ratio used to rescale all length with the right amount
            for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it) {
                target_length[*v_it]= target_length[*v_it]*ratio ; // rescale all length
            }
        }
    }

	void MeshProcessing::split_long_edges()
	{
		Mesh::Edge_iterator     e_it, e_end(mesh_.edges_end());
		Mesh::Vertex   v0, v1, v;
		bool            finished;
        int             i,j;

		Mesh::Vertex_property<Point> normals = mesh_.vertex_property<Point>("v:normal");
		Mesh::Vertex_property<Scalar> target_length = mesh_.vertex_property<Scalar>("v:length", 0);

        cout << "Splitting long edges..." << endl;//DEBUGGING

		for (finished = false, i = 0; !finished && i < 100; ++i)
		{
            j=1;
			finished = true;
			// ------------- IMPLEMENT HERE ---------
			// INSERT CODE:
			//  Compute the desired length as the mean between the property target_length of two vertices of the edge
			//  If the edge is longer than 4/3 * desired length
			//		add the midpoint to the mesh
			//		set the interpolated normal and interpolated vtargetlength_ property to the vertex
			//		split the edge with this vertex (use openMesh function split)
			// Leave the loop running until no splits are done (use the finished variable)
			// ------------- IMPLEMENT HERE ---------

            //Loop over all edges on the mesh
            Scalar target_l = 0; //for clarity
            Scalar l = 0; // for clarity
            for(e_it = mesh_.edges_begin(); e_it != e_end; ++e_it){

                //Taking the vertices of the edge e_it
                v0 = mesh_.to_vertex(mesh_.halfedge(*e_it,0));
                v1 = mesh_.to_vertex(mesh_.halfedge(*e_it,1));

                l = mesh_.edge_length(*e_it);
                target_l = 0.5 * (target_length[v0] + target_length[v1]);

                //Checking if edge length is greater than the 4/3 of the mean of the target length of each vertex
                if(l > (4/3) * target_l){

                    //Setting the new vertex replacing the two previous ones
                    v = mesh_.add_vertex((mesh_.position(v0) + mesh_.position(v1))/2);
                    normals[v] = (normals[v0] + normals[v1])/2;
                    target_length[v] = (target_length[v0] + target_length[v1])/2;

//                    cout << "Edge length : " << l << endl;//DEBUGGING
//                    cout << "Target length : " << dl << endl;//DEBUGGING
//                    cout << "Splitting the "<< j<< "th edge on iterations : " << i << endl; //DEBUGGING

                    //splitting with the vertex v
                    mesh_.split(*e_it,v);

                    //to continue the work, if no splitting finished will stay true
                    finished = false;
                j++; //DEBUGGING
                }
            }
            cout << j << " edges splitted on iteration " << i << endl; //DEBUGGING
		}

        cout << "==================================" << endl;//DEBUGGING
	}

    void MeshProcessing::collapse_short_edges()
    {
        Mesh::Edge_iterator     e_it, e_end(mesh_.edges_end());
        Mesh::Vertex   v0, v1;
        Mesh::Halfedge  h01, h10;
        bool            finished, b0, b1;
        int             i,j;
        bool            hcol01, hcol10;

        Mesh::Vertex_property<Scalar> target_length = mesh_.vertex_property<Scalar>("v:length", 0);

        cout << "Collapsing short edges..." << endl;//DEBUGGING


        for (finished = false, i = 0; !finished && i < 100; ++i)
        {
            finished = true;

            for (e_it = mesh_.edges_begin(); e_it != e_end; ++e_it)
            {
                if (!mesh_.is_deleted(*e_it)) // might already be deleted
                {
                    j = 1;
                    // ------------- IMPLEMENT HERE ---------
                    // INSERT CODE:
                    // Compute the desired length as the mean between the property vtargetlength_ of two vertices of the edge
                    // If the edge is shorter than 4/5 of the desired length
                    //		Check if halfedge connects a boundary vertex with a non-boundary vertex. If so, don't collapse.
                    //		Check if halfedges collapsible
                    //		Select the halfedge to be collapsed if at least one halfedge can be collapsed
                    //		Collapse the halfedge
                    // Leave the loop running until no collapse has been done (use the finished variable)
                    // ------------- IMPLEMENT HERE ---------

                    //Taking halfedges and vertices of the edge e_it
                    h01 = mesh_.halfedge(*e_it,0);
                    h10 = mesh_.halfedge(*e_it,1);
                    v0 = mesh_.to_vertex(h01);
                    v1 = mesh_.to_vertex(h10);

                    // Checking if edge length is greater than the 4/5 of the mean of the target length of each vertex
                    if(mesh_.edge_length(*e_it) < (2./5.) * (target_length[v0] + target_length[v1])){
                        finished = false;
                        if (mesh_.is_collapse_ok(h01) && mesh_.is_collapse_ok(h10) &&
                                !mesh_.is_boundary(v0) && !mesh_.is_boundary(v1)) {
                            if (mesh_.valence(v0) < mesh_.valence(v1)) {
                                mesh_.collapse(h10); ++j;
                            } else {
                                mesh_.collapse(h01); ++j;
                            }
                        } else if (mesh_.is_collapse_ok(h01)){
                            mesh_.collapse(h01); ++j;
                        } else if (mesh_.is_collapse_ok(h10)){
                            mesh_.collapse(h10); ++j;
                        } else {
                            finished = true;
                        }
                    }
                }
            }
           cout << j << " edges collapsed on iteration " << i << endl; //DEBUGGING
        }

        mesh_.garbage_collection();

        if (i == 100) std::cerr << "collapse break\n";
        cout << "==================================" << endl;//DEBUGGING
    }

	void MeshProcessing::equalize_valences()
	{
		Mesh::Edge_iterator     e_it, e_end(mesh_.edges_end());
		Mesh::Vertex   v0, v1, v2, v3;
		Mesh::Halfedge   h;
		int             val0, val1, val2, val3;
		int             val_opt0, val_opt1, val_opt2, val_opt3;
		int             ve0, ve1, ve2, ve3, ve_before, ve_after;
		bool            finished;
        int             i,j;

        cout << "Equalizing valences... " << endl;
		// flip all edges
		for (finished = false, i = 0; !finished && i < 100; ++i)
		{
			finished = true;
            j = 0; //DEBUGGING

			for (e_it = mesh_.edges_begin(); e_it != e_end; ++e_it)
			{
				if (!mesh_.is_boundary(*e_it))
				{
					// ------------- IMPLEMENT HERE ---------
					//  Extract valences of the four vertices involved to an eventual flip.
					//  Compute the sum of the squared valence deviances before flip
					//  Compute the sum of the squared valence deviances after an eventual flip
					//  If valence deviance is decreased and flip is possible, flip the vertex
					//  Leave the loop running until no collapse has been done (use the finished variable)
                    // ------------- IMPLEMENT HERE ---------

                  // defining the four vertices in the flip process
                  // The originally connected vertices are v0 and v1
                  h = mesh_.halfedge(*e_it,0);  // h0 connects v1 to v0
                  v0 = mesh_.to_vertex(h);      // v0
                  h = mesh_.next_halfedge(h);   // h2 connects v0 to v2
                  v2 = mesh_.to_vertex(h);      // v2

                  h = mesh_.halfedge(*e_it,1);  // h1 connects v0 to v1
                  v1 = mesh_.to_vertex(h);      // v1
                  h = mesh_.next_halfedge(h);   // h3 connects v1 to v3
                  v3 = mesh_.to_vertex(h);      // v3

                  //the actual valences
                  val0 = mesh_.valence(v0);
                  val1 = mesh_.valence(v1);
                  val2 = mesh_.valence(v2);
                  val3 = mesh_.valence(v3);

                  //the optimal valences (dependant of the boundary)
                  if (!mesh_.is_boundary(v0)){  val_opt0 = 6;   }
                  else                       {  val_opt0 = 4;   }

                  if (!mesh_.is_boundary(v1)){  val_opt1 = 6;   }
                  else                       {  val_opt1 = 4;   }

                  if (!mesh_.is_boundary(v2)){  val_opt2 = 6;   }
                  else                       {  val_opt2 = 4;   }

                  if (!mesh_.is_boundary(v3)){  val_opt3 = 6;   }
                  else                       {  val_opt3 = 4;   }

                  //the squared valence excess before the flip
                  ve0 = pow(val0 - val_opt0,2);
                  ve1 = pow(val1 - val_opt1,2);
                  ve2 = pow(val2 - val_opt2,2);
                  ve3 = pow(val3 - val_opt3,2);

                  //total squared valence excess before the flip
                  ve_before = ve0 + ve1 + ve2 + ve3;

                  //the squared valence excess after the flip
                  ve0 = pow(val0-1 - val_opt0,2);
                  ve1 = pow(val1-1 - val_opt1,2);
                  ve2 = pow(val2+1 - val_opt2,2);
                  ve3 = pow(val3+1 - val_opt3,2);

                  //total squared valence excess after the flip
                  ve_after = ve0 + ve1 + ve2 + ve3;

                  //if the valence excess decreases with the flip
                  //and the flip is possible we flip it
                  if( ve_after < ve_before && mesh_.is_flip_ok(*e_it))
                  {
                      mesh_.flip(*e_it);
                      finished = false;
                      ++j;
                  }
				}
			}
            cout << j << " edges flipped on iteration " << i << endl;
		}

		if (i == 100) std::cerr << "flip break\n";
     cout << "==================================" << endl;//DEBUGGING
	}

    void MeshProcessing::tangential_relaxation()
    {
        Mesh::Vertex_iterator     v_it, v_end(mesh_.vertices_end()); // Elle pose des problèmes donc je l'ai mise en commentaire
        Mesh::Vertex_around_vertex_circulator   vv_c, vv_end;
        int    valence;
        Point     u, n;
        Point     laplace;
        Scalar lambda = 0.01;

        Mesh::Vertex_property<Point> normals = mesh_.vertex_property<Point>("v:normal");
        Mesh::Vertex_property<Point> update = mesh_.vertex_property<Point>("v:update");

                cout << "Tangential smoothing ... " << endl;

        // smooth
        for (int iters = 0; iters < 10; ++iters)
        {
            for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it)
            {
                if (!mesh_.is_boundary(*v_it))
                {
                    // ------------- IMPLEMENT HERE ---------
                    //  Compute uniform laplacian curvature approximation vector
                    //  Compute the tangential component of the laplacian vector and move the vertex
                    //  Store smoothed vertex location in the update vertex property.
                    // ------------- IMPLEMENT HERE ---------
                    vv_c = mesh_.vertices(*v_it);
                    vv_end = vv_c;

                    valence = 0;
                    laplace = Point(0);

                    //computing uniform laplacian
                    do {
                        laplace += mesh_.position(*vv_c) - mesh_.position(*v_it);
                        ++valence;
                    } while(++vv_c != vv_end);
                    laplace /= valence;

                    //taking only the tangential part L_t = L - L_n and L_n = (L*n)*n
                    n = normals[*v_it];
                    u = laplace - dot(laplace,n)*n;
                    update[*v_it] = mesh_.position(*v_it) + lambda * u;
                }
            }
                    //updating the new positions
                     for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it)
                         if (!mesh_.is_boundary(*v_it))
                            mesh_.position(*v_it) = update[*v_it];
        }
        cout << "done" << endl;
        cout << "==================================" << endl;//DEBUGGING
        }


	void MeshProcessing::calc_uniform_mean_curvature() {
		Mesh::Vertex_property<Scalar> v_unicurvature =
			mesh_.vertex_property<Scalar>("v:unicurvature", 0.0f);
		// ------------- IMPLEMENT HERE ---------
		// For each non-boundary vertex, approximate mean curvature using
		// the length of the uniform Laplacian approximation
		// Save your approximation in unicurvature vertex property of the mesh.
		// ------------- IMPLEMENT HERE ---------
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

		// ------------- IMPLEMENT HERE ---------
		// For all non-boundary vertices, approximate the mean curvature using
		// the length of the Laplace-Beltrami approximation.
		// Save your approximation in v_curvature vertex property of the mesh.
		// Use the weights from calc_weights(): e_weight and v_weight
        // ------------- IMPLEMENT HERE ---------
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

		// ------------- IMPLEMENT HERE ---------
		// For each non-boundary vertex, approximate Gaussian curvature,
		// and store it in the vertex property v_gauss_curvature.
		// Hint: When calculating angles out of cross products make sure the value
		// you pass to the acos function is between -1.0 and 1.0.
		// Use the v_weight property for the area weight.
        // ------------- IMPLEMENT HERE ---------
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

	void MeshProcessing::calc_weights() {
		calc_edges_weights();
		calc_vertices_weights();
	}

/*	void MeshProcessing::calc_edges_weights() {
		auto e_weight = mesh_.edge_property<Scalar>("e:weight", 0.0f);
		auto points = mesh_.vertex_property<Point>("v:point");

		Mesh::Halfedge h0, h1, h2;
		Point p0, p1, p2, d0, d1;

		for (auto e : mesh_.edges())
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
				e_weight[e] += dot(d0, d1) / norm(cross(d0, d1));
			}

			if (!mesh_.is_boundary(h1))
			{
				h2 = mesh_.next_halfedge(h1);
				p2 = points[mesh_.to_vertex(h2)];
				d0 = p0 - p2;
				d1 = p1 - p2;
				e_weight[e] += dot(d0, d1) / norm(cross(d0, d1));
			}
		}
    }*/

	void MeshProcessing::calc_vertices_weights() {
		Mesh::Face_around_vertex_circulator vf_c, vf_end;
		Mesh::Vertex_around_face_circulator fv_c;
		Scalar area;
		auto v_weight = mesh_.vertex_property<Scalar>("v:weight", 0.0f);

		for (auto v : mesh_.vertices()) {
			area = 0.0;
			vf_c = mesh_.faces(v);

			if (!vf_c) {
				continue;
			}

			vf_end = vf_c;

			do {
				fv_c = mesh_.vertices(*vf_c);

				const Point& P = mesh_.position(*fv_c);  ++fv_c;
				const Point& Q = mesh_.position(*fv_c);  ++fv_c;
				const Point& R = mesh_.position(*fv_c);

				area += norm(cross(Q - P, R - P)) * 0.5f * 0.3333f;

			} while (++vf_c != vf_end);

			v_weight[v] = 0.5 / area;
		}
	}

	void MeshProcessing::load_mesh(const string &filename) {
		if (!mesh_.read(filename)) {
			std::cerr << "Mesh not found, exiting." << std::endl;
			exit(-1);
		}

		cout << "Mesh " << filename << " loaded." << endl;
		cout << "# of vertices : " << mesh_.n_vertices() << endl;
		cout << "# of faces : " << mesh_.n_faces() << endl;
		cout << "# of edges : " << mesh_.n_edges() << endl;

		// Compute the center of the mesh
		mesh_center_ = Point(0.0f, 0.0f, 0.0f);
		for (auto v : mesh_.vertices()) {
			mesh_center_ += mesh_.position(v);
		}
		mesh_center_ /= mesh_.n_vertices();

		// Compute the maximum distance from all points in the mesh and the center
		dist_max_ = 0.0f;
		for (auto v : mesh_.vertices()) {
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
		for (auto v : mesh_.vertices()) {
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
		color_coding(vertex_valence, &mesh_, v_color_valence, 3 /* min */,
			8 /* max */);
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

		for (auto f : mesh_.faces()) {
			std::vector<float> vv(3);
			int k = 0;
			for (auto v : mesh_.vertices(f)) {
				vv[k] = v.idx();
				++k;
			}
			indices_.col(j) << vv[0], vv[1], vv[2];
			++j;
		}

		j = 0;
		for (auto v : mesh_.vertices()) {
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
		Mesh::Vertex_property<Color> color_prop, Scalar min_value,
		Scalar max_value, int bound) {
		// Get the value array
		std::vector<Scalar> values = prop.vector();

		if (min_value == 0.0 && max_value == 0.0) {
			// discard upper and lower bound
			unsigned int n = values.size() - 1;
			unsigned int i = n / bound;
			std::sort(values.begin(), values.end());
			min_value = values[i];
			max_value = values[n - 1 - i];
		}

		// map values to colors
		for (auto v : mesh->vertices())
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
		v0 = min_value + 0.0 / 4.0 * (max_value - min_value);
		v1 = min_value + 1.0 / 4.0 * (max_value - min_value);
		v2 = min_value + 2.0 / 4.0 * (max_value - min_value);
		v3 = min_value + 3.0 / 4.0 * (max_value - min_value);
		v4 = min_value + 4.0 / 4.0 * (max_value - min_value);

		Color col(1.0f, 1.0f, 1.0f);

		if (value < v0) {
			col = Color(0, 0, 1);
		}
		else if (value > v4) {
			col = Color(1, 0, 0);
		}
		else if (value <= v2) {
			if (value <= v1) { // [v0, v1]
				Scalar u = (value - v0) / (v1 - v0);
				col = Color(0, u, 1);
			}
			else { // ]v1, v2]
				Scalar u = (value - v1) / (v2 - v1);
				col = Color(0, 1, 1 - u);
			}
		}
		else {
			if (value <= v3) { // ]v2, v3]
				Scalar u = (value - v2) / (v3 - v2);
				col = Color(u, 1, 0);
			}
			else { // ]v3, v4]
				Scalar u = (value - v3) / (v4 - v3);
				col = Color(1, 1 - u, 0);
			}
		}
		return col;
	}


}


