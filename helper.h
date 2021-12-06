#pragma once
#include <vector>
#include <set>
#include <map>
#include "TinyVector.h"
#include "Tree.h"

void greedy_graph_coloring(const size_t num_vertices, const std::vector<std::set<size_t>>& edges, std::vector<std::vector<size_t>>& colored_vertices);

void get_graph_component(const std::vector<std::set<size_t>>& edges, std::vector<std::set<size_t>>& components);

template<typename T>
T compute_tri_area(const TinyVector<T, 3>& v0, const TinyVector<T, 3>& v1, const TinyVector<T, 3>& v2)
{
	return ((v1 - v0).Cross(v2 - v0)).Length() * (T)0.5;
}

void compute_vert_mean_curvature(const std::vector<TinyVector<double, 3>>& pos, const std::vector<TinyVector<size_t, 3>> &faces, std::vector<double>& curvature);

void update_tree_color(const std::vector<size_t>& new_color, TreeNode<size_t>* t);

bool get_tree_from_convex_graph(const std::vector<std::set<size_t>> &graph, const std::map<std::pair<size_t, size_t>, int>& flag_fpconvex, bool flag_convex, TreeNode<size_t>* tn, int layer);

TinyVector<double, 3> perturb_normal(const TinyVector<double, 3> normal, double angle_noise1, double angle_noise2);

bool sample_pts_from_mesh(const std::vector<TinyVector<double, 3>>& tri_verts, const std::vector<TinyVector<size_t, 3>>& tri_faces, const std::vector<TinyVector<double, 3>>& tri_normals, const std::vector<int>& tri_face_masks, int n_sample, std::vector<TinyVector<double, 3>>& output_pts, std::vector<TinyVector<double, 3>>& output_normals, std::vector<int>& output_masks, double sigma = -1.0, double sigma_n = -1.0);

#include "Mesh3D.h"
//cluster mesh faces according to feature edges
void get_grouped_edges(Mesh3d& m, const std::vector<bool>& he_feautre_flag, const std::vector<std::vector<int>>& feature_v2he, std::vector<std::vector<int>>& grouped_features, std::vector<int>& he2gid);
int cluster_mesh_faces(Mesh3d* m, const std::vector<bool>& he_feature_flag, std::vector<std::vector<int>> &grouped_features, int cluster_begin_id, std::vector<int>& face2cluster, std::vector<std::pair<int, int>>& feature_twoface_colors);

void get_cluster_from_coloring(const std::vector<int>& face_color, int color_start_id, std::vector<std::vector<int>>& face_clusters);
int merge_clusters(Mesh3d* m, const std::vector<bool>& he_feature_flag, int cluster_start_id, int n_cluster, const std::vector<std::pair<int, int>>& feature_twoface_colors, std::vector<int>& face_color);