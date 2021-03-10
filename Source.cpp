﻿#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <random>
#include <algorithm>
#include <direct.h>
#include "cxxopts.hpp"
#include "Mesh3D.h"
#include "helper.h"
#include "happly.h"
#include "Tree.h"

#define DENOMINATOR_EPS 1e-6

using namespace MeshLib;

std::string GetFileExtension(const std::string& FileName)
{
	if (FileName.find_last_of(".") != std::string::npos)
		return FileName.substr(FileName.find_last_of(".") + 1);
	return "";
}

bool load_feature_file(const char* filename, std::vector<std::pair<int, int>>& ungrouped_feature)
{
	ungrouped_feature.clear();
	std::ifstream ifs(filename);
	if (!ifs.is_open())
	{
		std::cout << "Cannot Open Input Feature File" << std::endl;
		return false;
	}

	//fea file: first part
	int pair_size = 0;
	ifs >> pair_size;
	std::pair<int, int> tmp_pair(-1, -1);
	for (size_t i = 0; i < pair_size; i++)
	{
		ifs >> tmp_pair.first >> tmp_pair.second;
		ungrouped_feature.push_back(tmp_pair);
	}

	//not consider grouped feature
	return true;

}

void save_feature_file(const char* filename, const std::vector<std::pair<int, int>>& ungrouped_feature)
{
	std::ofstream ofs(filename);
	ofs << ungrouped_feature.size() << std::endl;
	for (size_t i = 0; i < ungrouped_feature.size(); i++)
	{
		ofs << ungrouped_feature[i].first << " " << ungrouped_feature[i].second << std::endl;
	}

	ofs.close();
}

void save_conf_file(const char* filename, const std::string str, bool flag_convex = true)
{
	std::ofstream ofs(filename);
	ofs << "csg{\n    list = ";
	ofs << str << std::endl;
	ofs << "    flag_convex = " << int(flag_convex) << "," << std::endl;
	ofs << "}";
	
	ofs.close();
}

bool check_mesh_edge_convex(Mesh3d* m, HE_edge<double>* he)
{
	int hetri_vertidsum = 0, hepairtri_vertidsum = 0;
	
	HE_edge<double>* begin_edge = he;
	HE_edge<double>* edge = he;
	do
	{
		hetri_vertidsum += edge->vert->id;
		edge = edge->next;
	} while (edge != begin_edge);

	begin_edge = he->pair;
	edge = he->pair;
	do
	{
		hepairtri_vertidsum += edge->vert->id;
		edge = edge->next;
	} while (edge != begin_edge);

	int hetri_otherid = hetri_vertidsum - he->vert->id - he->pair->vert->id;
	int hepairtri_otherid = hepairtri_vertidsum - he->vert->id - he->pair->vert->id;
	double product = he->face->normal.Dot(m->get_vertices_list()->at(hetri_otherid)->pos - m->get_vertices_list()->at(hepairtri_otherid)->pos);
	if (product > 0.0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool dijkstra_mesh(Mesh3d* m, size_t start_vert, std::vector<size_t>& prev_map, std::vector<double>& dist)
{
	prev_map.clear();
	dist.clear();
	dist.resize(m->get_num_of_vertices(), DBL_MAX);
	prev_map.resize(m->get_num_of_vertices(), 0);
	dist[start_vert] = 0.0;
	std::set<std::pair<double, HE_vert<double>*>> vqueue;
	for (size_t i = 0; i < m->get_num_of_vertices(); i++)
	{
		HE_vert<double>* hv = m->get_vertex(i);
		vqueue.insert(std::pair<double, HE_vert<double>*>(dist[hv->id], hv));
	}
	
	while (!vqueue.empty())
	{
		std::set<std::pair<double, HE_vert<double>*>>::iterator iter = vqueue.begin();
		HE_vert<double>* u = iter->second;
		vqueue.erase(iter);
		if (u->edge == NULL)
			continue;

		HE_edge<double>* he = u->edge;
		do
		{
			double alt = dist[u->id] + (u->pos - he->vert->pos).Length();
			if (alt < dist[he->vert->id])
			{
				vqueue.erase(std::pair<double, HE_vert<double>*>(dist[he->vert->id], he->vert));
				dist[he->vert->id] = alt;
				prev_map[he->vert->id] = u->id;
				vqueue.insert(std::pair<double, HE_vert<double>*>(dist[he->vert->id], he->vert));
			}
			he = he->pair->next;
		} while (he != u->edge);
	}

	return true;
}

TinyVector<double, 3> perturb_normal(const TinyVector<double, 3> normal, double angle_noise1, double angle_noise2)
{
	//get angle from normal
	//assume normal is normalized
	double theta = acos(normal[2]);
	double phi = acos(normal[0] / (sqrt(1 - normal[2] * normal[2]) + DENOMINATOR_EPS));
	double phi_ref = asin(normal[1] / (sqrt(1 - normal[2] * normal[2]) + DENOMINATOR_EPS));
	if (phi_ref < 0)
		phi = 2 * M_PI - phi;
	theta = theta + angle_noise1 / 180.0 * M_PI;
	phi = phi + angle_noise2 / 180.0 * M_PI;
	TinyVector<double, 3> normal_new;
	normal_new[0] = sin(theta) * cos(phi);
	normal_new[1] = sin(theta) * sin(phi);
	normal_new[2] = cos(theta);
	return normal_new;
}


int main(int argc, char** argv)
{
	
	//select model firstly: 0 for normalization and 1 for sampling
	try
	{
		cxxopts::Options options("FeaturedModelPointSample", "Point Sampling program for featured CAD models (author: Haoxiang Guo, Email: guohaoxiangxiang@gmail.com)");
		options
			.positional_help("[optional args]")
			.show_positional_help()
			.allow_unrecognised_options()
			.add_options()
			("i,input", "input mesh (obj/off format)", cxxopts::value<std::string>())
			("f,feature", "input feature file (fea format)", cxxopts::value<std::string>())
			("o,output", "output mesh/points (obj/off/points/xyz format)", cxxopts::value<std::string>())
			("k,mask", "output mask file (txt format)", cxxopts::value<std::string>())
			("fs", "number of samples on feature edges(default: 10000)", cxxopts::value<int>())
			("ns", "number of samples on non-feature faces(default: 40000)", cxxopts::value<int>())
			("m,mode", "processing mode: 0 for normalization and 1 for feature sample", cxxopts::value<int>())
			("c,color", "whether coloring is used, 0: not used, 1: used, default: 0", cxxopts::value<int>())
			("mp", "maximum number of patches in each colored cluster, only work for csg, default -1(no upper bound)", cxxopts::value<int>())
			("cot", "whether cotangent weight is used for sampling, 0: not used, 1: used, default: 0", cxxopts::value<int>())
			("s,sigma", "sigma for noisy points position, default 0.0", cxxopts::value<double>())
			("sn", "sigma for noisy points normal in degrees, default 0.0", cxxopts::value<double>())
			("csg", "whether generating csg tree for model, default: 0", cxxopts::value<int>())
			("convex", "whether the first layer is convex, default: 0", cxxopts::value<int>())
			("r,repair", "whether the feature edges are repaired, default: 1", cxxopts::value<int>())
			("h,help", "print help");
		
		auto result = options.parse(argc, argv);
		if (result.count("help"))
		{
			std::cout << options.help({ "", "Group" }) << std::endl;
			exit(0);
		}
		int n_nonfeature_sample = 40000;
		int n_feature_sample = 10000;
		double sigma = 0.0;
		double sigma_n = 0.0;
		assert(result.count("m"));
		int processing_mode = result["m"].as<int>();
		assert(result.count("i") && result.count("o"));
		auto& inputfile = result["i"].as<std::string>();
		auto& outputfile = result["o"].as<std::string>();
		//output pts by colors
		int last_dot = (int)outputfile.find_last_of(".");
		auto output_prefix = outputfile.substr(0, last_dot);
		int flag_csg = 0;

		std::string inputext = GetFileExtension(inputfile);
		Mesh3d mesh;
		if (inputext == "obj")
			mesh.load_obj(inputfile.c_str());
		else if (inputext == "off")
			mesh.load_off(inputfile.c_str());
		std::cout << "verts: " << mesh.get_vertices_list()->size() << " face:  " << mesh.get_faces_list()->size() << std::endl;

		if (processing_mode == 0)
		{
			//normalization part begin
			//[-0.9, 9]^3
			std::vector<TinyVector<double, 3>> pts_nl(mesh.get_vertices_list()->size());
			double max_range = mesh.xmax - mesh.xmin;
			max_range = max_range < (mesh.ymax - mesh.ymin) ? (mesh.ymax - mesh.ymin) : max_range;
			max_range = max_range < (mesh.zmax - mesh.zmin) ? (mesh.zmax - mesh.zmin) : max_range;

			double xcenter = (mesh.xmin + mesh.xmax) / 2;
			double ycenter = (mesh.ymin + mesh.ymax) / 2;
			double zcenter = (mesh.zmin + mesh.zmax) / 2;
			std::cout << "center " << xcenter << " " << ycenter << " " << zcenter << std::endl;

			for (size_t i = 0; i < mesh.get_vertices_list()->size(); i++)
			{
				mesh.get_vertices_list()->at(i)->pos[0] = (mesh.get_vertices_list()->at(i)->pos[0] - xcenter) / max_range * 1.8;
				mesh.get_vertices_list()->at(i)->pos[1] = (mesh.get_vertices_list()->at(i)->pos[1] - ycenter) / max_range * 1.8;
				mesh.get_vertices_list()->at(i)->pos[2] = (mesh.get_vertices_list()->at(i)->pos[2] - zcenter) / max_range * 1.8;
			}

			//output mesh
			std::string outputext = GetFileExtension(outputfile);
			if (outputext == "obj")
			{
				mesh.write_obj(outputfile.c_str());

			}
			else if (outputext == "off")
			{
				mesh.write_off(outputfile.c_str());
			}
			return 1;
		}
		else if (processing_mode == 1)
		{
			//first sample feature parts then non-feature parts
			//mask: 
			//feature: 0
			//non feature: 1,2,3...indicating coloring
			assert(result.count("f") && result.count("k"));
			auto& inputfeaturefile = result["f"].as<std::string>();
			auto& outputmaskfile = result["k"].as<std::string>();
			if (result.count("fs"))
				n_feature_sample = result["fs"].as<int>();
			if (result.count("ns"))
				n_nonfeature_sample = result["ns"].as<int>();
			if (result.count("s"))
				sigma = result["s"].as<double>();
			if (result.count("sn"))
				sigma_n = result["sn"].as<double>();
			bool flag_repair_features = true;
			if (result.count("r"))
				flag_repair_features = result["r"].as<int>();

			//at least sample points with the same shape as face
			if (n_nonfeature_sample < mesh.get_num_of_faces())
			{
				n_nonfeature_sample = mesh.get_num_of_faces();
			}

			std::vector<int> sample_mask(n_feature_sample + n_nonfeature_sample, 0);
			
			std::vector<std::pair<int, int>> ungrouped_features;
			load_feature_file(inputfeaturefile.c_str(), ungrouped_features);
			std::vector<TinyVector<double, 3>> sample_pts, sample_pt_normals;

			//feature check: no hanging feature
			//std::vector<size_t> feature_degree_v(mesh.get_num_of_vertices(), 0);
			std::vector<std::vector<int>> feature_v2he(mesh.get_num_of_vertices()); //id of hes ematating from each vertex

			std::vector<std::pair<int, int>> ungrouped_features_new;
			for (size_t i = 0; i < ungrouped_features.size(); i++)
			{
				int id0 = ungrouped_features[i].first;
				int id1 = ungrouped_features[i].second;
				//feature_degree_v[id0]++;
				//feature_degree_v[id1]++;

				HE_edge<double>* begin_edge = mesh.get_vertices_list()->at(id0)->edge;
				HE_edge<double>* edge = mesh.get_vertices_list()->at(id0)->edge;
				bool flag_found = false;
				do
				{
					if (id1 == edge->vert->id)
					{
						feature_v2he[id0].push_back(edge->id);
						feature_v2he[id1].push_back(edge->pair->id);
						flag_found = true;
						break;
					}
					edge = edge->pair->next;
				} while (edge != begin_edge);
				//assert(flag_found == true);
				if (flag_found == true)
				{
					ungrouped_features_new.push_back(ungrouped_features[i]);
				}
			}

			if (ungrouped_features.size() != ungrouped_features_new.size())
				ungrouped_features = ungrouped_features_new;


			std::vector<size_t> turn_verts, hanging_verts;

			for (size_t i = 0; i < mesh.get_num_of_vertices(); i++)
			{
				//assert(feature_degree_v[i] == feature_v2he[i].size());
				if (feature_v2he[i].size() == 1)
				{
					std::cout << "input file: " << inputfile << std::endl;
					std::cout << "hanging vertex exists: " << i << std::endl;
					//turn_verts.push_back(i);
					hanging_verts.push_back(i);
					//return 0;
				}
				else if (feature_v2he[i].size() == 2)
				{
					//check edge convex status
					bool flag_convex_edge0 = check_mesh_edge_convex(&mesh, mesh.get_edges_list()->at(feature_v2he[i][0]));
					bool flag_convex_edge1 = check_mesh_edge_convex(&mesh, mesh.get_edges_list()->at(feature_v2he[i][1]));
					if (flag_convex_edge0 != flag_convex_edge1)
					{
						std::cout << "input file: " << inputfile << std::endl;
						std::cout << "turn vertex exists: " << i << std::endl;
						turn_verts.push_back(i);
					}
				}
			}

			//do not handle hanging vertex
			if (!hanging_verts.empty())
			{
				return 0;
			}

			/*if (!flag_repair_features && !(turn_verts.size() + hanging_verts.size() == 0))
			{
				return 0;
			}*/

			//feature parts first
			std::random_device rd;
			std::mt19937 e2(rd());
			std::uniform_real_distribution<double> unif_dist(0, 1);
			//std::normal_distribution<double> normal_dist(0, sigma);
			std::uniform_real_distribution<double> normal_dist(-sigma, sigma);
			std::uniform_real_distribution<double> angle_unif_dist(-sigma_n, sigma_n);

			std::vector<double> feature_length(ungrouped_features.size(), 0.0);
			double total_feature_length = 0.0;
			for (size_t i = 0; i < ungrouped_features.size(); i++)
			{
				int id0 = ungrouped_features[i].first;
				int id1 = ungrouped_features[i].second;
				feature_length[i] = (mesh.get_vertices_list()->at(id0)->pos - mesh.get_vertices_list()->at(id1)->pos).Length();
				total_feature_length += feature_length[i];
			}
			
			std::vector<double> line_bound(ungrouped_features.size() + 1, 0.0);
			for (size_t i = 0; i < ungrouped_features.size(); i++)
			{
				line_bound[i + 1] = line_bound[i] + feature_length[i] / total_feature_length;
			}

			//sampling
			for (size_t i = 0; i < n_feature_sample; i++)
			{
				double u = unif_dist(e2);
				auto iter = std::upper_bound(line_bound.begin(), line_bound.end(), u);
				int fid = (int)std::distance(line_bound.begin(), iter);
				assert(fid != ungrouped_features.size() + 1);
				fid = std::max(0, fid - 1);
				//sample
				int id0 = ungrouped_features[fid].first;
				int id1 = ungrouped_features[fid].second;
				double s = unif_dist(e2);
				sample_pts.push_back(s * mesh.get_vertices_list()->at(id0)->pos + (1.0 - s) * mesh.get_vertices_list()->at(id1)->pos);
				sample_pt_normals.push_back(TinyVector<double, 3>(1.0, 0.0, 0.0));
			}

			std::vector<double> tri_area(mesh.get_faces_list()->size(), 0.0);
			double total_tri_area = 0.0;
			std::vector<TinyVector<size_t, 3>> tri_verts(mesh.get_faces_list()->size());
			//get tri list
			for (size_t i = 0; i < mesh.get_faces_list()->size(); i++)
			{
				HE_edge<double>* begin_edge = mesh.get_faces_list()->at(i)->edge;
				HE_edge<double>* edge = mesh.get_faces_list()->at(i)->edge;
				int local_id = 0;
				do
				{
					tri_verts[i][local_id++] = edge->pair->vert->id;
					edge = edge->next;
				} while (edge != begin_edge);
			}
			std::vector<TinyVector<double, 3>> vert_pos;
			for (size_t i = 0; i < mesh.get_vertices_list()->size(); i++)
			{
				vert_pos.push_back(mesh.get_vertices_list()->at(i)->pos);
			}

			//cluster faces: assuming the features are all close
			std::vector<bool> he_feature_flag(mesh.get_edges_list()->size(), false);
			for (size_t i = 0; i < ungrouped_features.size(); i++)
			{
				int id0 = ungrouped_features[i].first;
				int id1 = ungrouped_features[i].second;
				//iterate over all verts emanating from id0
				HE_edge<double>* edge = mesh.get_vertices_list()->at(id0)->edge;
				do
				{
					if (edge->vert->id == id1)
					{
						break;
					}
					edge = edge->pair->next;
				} while (edge != mesh.get_vertices_list()->at(id0)->edge);
				assert(edge->vert->id == id1);
				he_feature_flag[edge->id] = true;
				he_feature_flag[edge->pair->id] = true;
			}
			
			//repairing conducted here, he_feature_flag, ungrouped_features should be updated
			
			std::vector<std::vector<int>> grouped_features; //grouped features: only one he of a pair is stored
			std::vector<int> he2gid;
			get_grouped_edges(mesh, he_feature_flag, feature_v2he, grouped_features, he2gid);
			if (flag_repair_features && !turn_verts.empty())
			//if (true)
			{
				//group features first
				
				std::cout << "feature group size: " << grouped_features.size() << std::endl;
				//repair features
				std::vector<bool> flag_feature_points(mesh.get_vertices_list()->size(), false);
				for (size_t i = 0; i < he_feature_flag.size(); i++)
				{
					if (he_feature_flag[i])
					{
						flag_feature_points[mesh.get_edges_list()->at(i)->vert->id] = true;
					}
				}

				//repair turn vertex
				for (size_t i = 0; i < turn_verts.size(); i++)
				{
					size_t cur_vert = turn_verts[i];
					std::vector<bool> flag_feature_points_tmp = flag_feature_points;
					size_t cur_group = he2gid[feature_v2he[cur_vert][0]];
					for (auto heid : grouped_features[cur_group])
					{
						flag_feature_points_tmp[mesh.get_edges_list()->at(heid)->vert->id] = false;
						flag_feature_points_tmp[mesh.get_edges_list()->at(heid)->pair->vert->id] = false;
					}

					//distance from cur_vert to all other verts
					std::vector<size_t> prev_map;
					std::vector<double> dist;
					dijkstra_mesh(&mesh, cur_vert, prev_map, dist);
					
					/*size_t nnvid = 0;
					double shortest_dist = DBL_MAX;*/
					std::set<std::pair<double, size_t>> dist_id_set;
					for (size_t j = 0; j < flag_feature_points_tmp.size(); j++)
					{
						if (flag_feature_points_tmp[j])
						{
							/*if (dist[j] < shortest_dist)
							{
								shortest_dist = dist[j];
								nnvid = j;
							}*/
							dist_id_set.insert(std::pair<double, size_t>(dist[j], j));
						}
					}

					for (auto& dist_id : dist_id_set)
					{
						size_t nnvid = dist_id.second;
						bool flag_usable = true;
						std::vector<HE_edge<double>*> tmp_hes;
						//add path from nnvid to cur_vert
						while (nnvid != cur_vert)
						{
							size_t prev_vert = prev_map[nnvid];
							//ungrouped_features.push_back(std::pair<int, int>((int)nnvid, (int)prev_vert));
							HE_edge<double>* begin_edge = mesh.get_vertices_list()->at(nnvid)->edge;
							HE_edge<double>* edge = begin_edge;
							do
							{
								if (edge->vert->id == prev_vert)
								{
									break;
								}
								edge = edge->pair->next;
							} while (edge != begin_edge);
							assert(edge->vert->id == prev_vert);
							/*he_feature_flag[edge->id] = true;
							he_feature_flag[edge->pair->id] = true;*/
							if (he_feature_flag[edge->id])
							{
								flag_usable = false;
								break;
							}
							else
							{
								tmp_hes.push_back(edge);
							}

							nnvid = prev_vert;
						}
						if (flag_usable)
						{
							for (auto he : tmp_hes)
							{
								ungrouped_features.push_back(std::pair<int, int>(he->vert->id, he->pair->vert->id));
								he_feature_flag[he->id] = true;
								he_feature_flag[he->pair->id] = true;
								feature_v2he[he->vert->id].push_back(he->pair->id);
								feature_v2he[he->pair->vert->id].push_back(he->id);
							}
							break;
						}
					}

				
				}

				//save repaired features
				save_feature_file((output_prefix + "_repairturn.fea").c_str(), ungrouped_features);
			}
			
			std::vector<int> face_color_init(mesh.get_faces_list()->size(), -1); //starting from 1
			std::vector<std::vector<int>> face_clusters;
			std::vector<std::pair<int, int>> feature_twoface_colors;//color of faces on two sides of the feature, *.first < *.second


			//new version
			//int cluster_start_id = 1;
			////set face_color_init by propagation
			//get_grouped_edges(mesh, he_feature_flag, feature_v2he, grouped_features, he2gid);
			//int cluster_id = cluster_mesh_faces(&mesh, he_feature_flag, grouped_features, cluster_start_id, face_color_init, feature_twoface_colors);
			//std::cout << "cluster number before merging: " << cluster_id - 1 << std::endl;

			////get face_clusters
			////get_cluster_from_coloring(face_color_init, cluster_start_id, face_clusters);

			////merging: udpate face_color_init, cluster_id, face_clusters
			//cluster_id = merge_clusters(&mesh, he_feature_flag, cluster_start_id, cluster_id - 1,feature_twoface_colors, face_color_init);
			//std::cout << "cluster number after merging: " << cluster_id - 1 << std::endl;

			//new version end

#if 1
			//traditional version
			int start_id = 0;
			int cluster_id = 1;
			while (start_id != -1)
			{
				std::vector<int> onecluster;
				std::queue<int> q;
				q.push(start_id);
				face_color_init[start_id] = cluster_id;
				while (!q.empty())
				{
					int front = q.front();
					q.pop();
					onecluster.push_back(front);
					HE_edge<double>* edge = mesh.get_faces_list()->at(front)->edge;
					do
					{
						if (he_feature_flag[edge->id] == false)
						{
							//not feature
							int pair_fid = edge->pair->face->id;
							if (face_color_init[pair_fid] == -1)
							{
								q.push(pair_fid);
								face_color_init[pair_fid] = cluster_id;
							}
						}

						edge = edge->next;
					} while (edge != mesh.get_faces_list()->at(front)->edge);

				}

				face_clusters.push_back(onecluster);
				start_id = -1;
				//find next start_id
				for (size_t i = 0; i < face_color_init.size(); i++)
				{
					if (face_color_init[i] == -1)
					{
						start_id = i;
						break;
					}
				}
				cluster_id++;
			}
#endif 
			std::vector<int> face_color = face_color_init; //starting from 1
			int n_color = cluster_id;
			//coloring
			bool flag_coloring = false;
			if (result.count("c"))
				flag_coloring = result["c"].as<int>();

			if (result.count("csg"))
				flag_csg = result["csg"].as<int>();

			bool flag_first_convex = false;
			if (result.count("convex"))
				flag_first_convex = result["convex"].as<int>();

			if (flag_coloring && !flag_csg)
			{
				
				//need to be changed if csg is set as true
				std::vector<std::set<size_t>> connectivity(cluster_id - 1);
				//color - 1
				for (size_t i = 0; i < he_feature_flag.size(); i++)
				{
					if (he_feature_flag[i])
					{
						HE_edge<double>* e1 = mesh.get_edges_list()->at(i);
						HE_edge<double>* e2 = e1->pair;
						connectivity[face_color_init[e1->face->id] - 1].insert(face_color_init[e2->face->id] - 1);
					}
				}

				//print graph
				std::cout << "graph:" << std::endl;
				for (size_t i = 0; i < connectivity.size(); i++)
				{
					std::cout << i + 1 << ": ";
					for (auto v : connectivity[i])
					{
						std::cout << v + 1 << " ";
					}
					std::cout << std::endl;
				}

				std::vector<std::vector<size_t>> colored_vertices;
				greedy_graph_coloring(cluster_id - 1, connectivity, colored_vertices);
				std::cout << "number of colors: " << colored_vertices.size() << std::endl;
				n_color = colored_vertices.size() + 1;
				//update face_color
				for (size_t i = 0; i < colored_vertices.size(); i++)
				{
					for (size_t j = 0; j < colored_vertices[i].size(); j++)
					{
						size_t local_id = colored_vertices[i][j];
						for (size_t k = 0; k < face_clusters[local_id].size(); k++)
						{
							face_color[face_clusters[local_id][k]] = i + 1;
						}
					}
				}
			}

			if (flag_csg)
			{
				//coloring is based on vertices
				std::vector<std::set<size_t>> connectivity(cluster_id - 1);
				std::vector<std::set<size_t>> connectivity_v(cluster_id - 1); //connectivity based on vertices
				std::map<std::pair<size_t, size_t>, double> fp2product;
				std::map<std::pair<size_t, size_t>, bool> flag_fpconvex;
				//color - 1
				for (size_t i = 0; i < he_feature_flag.size(); i++)
				{
					if (he_feature_flag[i])
					{
						HE_edge<double>* e1 = mesh.get_edges_list()->at(i);
						HE_edge<double>* e2 = e1->pair;
						connectivity[face_color_init[e1->face->id] - 1].insert(face_color_init[e2->face->id] - 1);
						size_t fid1 = face_color_init[e1->face->id], fid2 = face_color_init[e2->face->id]; //starting from zero
						size_t minfid = std::min(fid1, fid2);
						size_t maxfid = std::max(fid1, fid2);
						
						//triangle face
						size_t tfid1 = e1->face->id, tfid2 = e2->face->id;
						size_t ev1 = e1->vert->id, ev2 = e2->vert->id;
						size_t tv1 = tri_verts[tfid1][0] + tri_verts[tfid1][1] + tri_verts[tfid1][2] - ev1 - ev2;
						size_t tv2 = tri_verts[tfid2][0] + tri_verts[tfid2][1] + tri_verts[tfid2][2] - ev1 - ev2;
						double product = e1->face->normal.Dot(vert_pos[tv2] - vert_pos[tv1]);

						std::pair<size_t, size_t> tmp_pair(minfid - 1, maxfid - 1);
						auto it = fp2product.find(tmp_pair);
						if (it == fp2product.end())
						{
							fp2product[tmp_pair] = product;
						}
						else
						{
							fp2product[tmp_pair] += product;
						}
					}
				}

				//init connectivity_v
				connectivity_v = connectivity;
				for (size_t i = 0; i < feature_v2he.size(); i++)
				{
					if (feature_v2he[i].size() > 2)
					{
						//vertex with degree larger than 3
						HE_edge<double>* ve_begin = mesh.get_vertices_list()->at(i)->edge;
						assert(ve_begin->pair->vert->id == i);
						HE_edge<double>* ve_iter = ve_begin;
						std::set<size_t> surounding_cs;
						do
						{
							/*int next_id = ve_iter->vert->id;
							if (next_id == id1)
							{
								te = ve_iter;
								break;
							}
							ve_iter = ve_iter->pair->next;*/
							surounding_cs.insert(face_color_init[ve_iter->face->id] -1);
							ve_iter = ve_iter->pair->next;
						} while (ve_iter != ve_begin);
						if (surounding_cs.size() != feature_v2he[i].size())
						{
							std::cout << "input file: " << inputfile << std::endl;
							std::cout << "feature wrong near: " << i << std::endl;
						}
						assert(surounding_cs.size() == feature_v2he[i].size());
						std::vector<size_t> surounding_cs_vector(surounding_cs.begin(), surounding_cs.end());
						for (size_t it = 0; it < surounding_cs_vector.size() - 1; it++)
						{
							size_t cur_fid = surounding_cs_vector[it];
							for (size_t it1 = it + 1; it1 < surounding_cs_vector.size(); it1++)
							{
								size_t n_fid = surounding_cs_vector[it1];
								connectivity_v[cur_fid].insert(n_fid);
								connectivity_v[n_fid].insert(cur_fid);
							}
						}
					}
				}
				

				for (auto& p : fp2product)
				{
					if (p.second < 0.0)
					{
						flag_fpconvex[p.first] = true;
					}
					else
					{
						flag_fpconvex[p.first] = false;
					}
				}

				//print graph
				/*std::cout << "dual graph ori:" << std::endl;
				for (size_t i = 0; i < connectivity.size(); i++)
				{
					std::cout << i << ": ";
					for (auto v : connectivity[i])
					{
						std::cout << v << " ";
					}
					std::cout << std::endl;
				}
				
				std::cout << "edge convex flag: " << std::endl;
				for (auto& p : flag_fpconvex)
				{
					std::cout << "edge: " << p.first.first << "-" << p.first.second << " : " << p.second << std::endl;
				}*/

				TreeNode<size_t> *tree = new TreeNode<size_t>;
				;
				bool flag_convex = true;
				//convex flag is set to true by default, but if the model contains multiple component, then it is set to concave

				std::vector<std::set<size_t>> components;
				get_graph_component(connectivity, components);
				if (components.size() == 1)
				{ 
					get_tree_from_convex_graph(connectivity, flag_fpconvex, true, tree, 0);
				}
				else
				{
					flag_convex = false;
					for (size_t ic = 0; ic < components.size(); ic++)
					{
						std::vector<std::set<size_t>> subgraph(connectivity.size());
						for (auto v : components[ic])
						{
							for (auto vn : connectivity[v])
							{
								if (std::find(components[ic].begin(), components[ic].end(), vn) != components[ic].end())
								{
									subgraph[v].insert(vn);
								}
							}
						}
						TreeNode<size_t>* child = new TreeNode<size_t>;
						get_tree_from_convex_graph(subgraph, flag_fpconvex, !flag_convex, child, 1);
						tree->children.push_back(child);
					}
				}

				//TreeNode<size_t>* tree_concave = new TreeNode<size_t>;
				//get_tree_from_convex_graph(connectivity, flag_fpconvex, false, tree_concave);
				////compare number of keys in the first layer
				//if (tree_concave->keys.size() > tree->keys.size() && !flag_first_convex)
				//{
				//	tree = tree_concave;
				//	flag_convex = false;
				//}
				
				std::cout << "convex status: " << flag_convex << std::endl;
				
				if (flag_coloring)
				{
					int max_patch_per_cluster = -1;
					if (result.count("mp"))
						max_patch_per_cluster = result["mp"].as<int>();
					std::vector<size_t> cluster_color(connectivity_v.size(), -1);
					n_color = tree_coloring<size_t>(tree, connectivity_v, cluster_color, 0, max_patch_per_cluster) + 1;
					//only color the vertex of root
					//std::vector<size_t> o2n(connectivity.size(), size_t(-1));
					//std::vector<size_t> n2o;
					//size_t ncount = 0;
					//for (auto v : tree->keys)
					//{
					//	o2n[v] = ncount++;
					//	n2o.push_back(v);
					//}
					//std::vector<std::set<size_t>> local_con(ncount);
					//for (auto v : tree->keys)
					//{
					//	for (auto vn : connectivity_v[v])
					//	{
					//		if (std::find(tree->keys.begin(), tree->keys.end(), vn) != tree->keys.end())
					//		{
					//			//both are found in keys
					//			local_con[o2n[v]].insert(o2n[vn]);
					//		}
					//	}
					//}
					//std::vector<std::vector<size_t>> colored_vertices;
					//greedy_graph_coloring(ncount, local_con, colored_vertices);
					////update colored_vertices
					//std::vector<size_t> cluster_color(connectivity.size()); //starting from 0
					//for (size_t i = 0; i < colored_vertices.size(); i++)
					//{
					//	for (size_t j = 0; j < colored_vertices[i].size(); j++)
					//	{
					//		colored_vertices[i][j] = n2o[colored_vertices[i][j]];
					//		cluster_color[colored_vertices[i][j]] = i;
					//	}
					//}
					//std::cout << "num of colors pre :" << colored_vertices.size() << std::endl;

					//for (size_t i = 0; i < connectivity.size(); i++)
					//{
					//	//append
					//	if (std::find(tree->keys.begin(), tree->keys.end(), i) == tree->keys.end())
					//	{
					//		cluster_color[i] = colored_vertices.size();
					//		std::vector<size_t> leaf_vert;
					//		leaf_vert.push_back(i);
					//		colored_vertices.push_back(leaf_vert);
					//	}
					//}
					//std::cout << "num of colors:" << colored_vertices.size() << std::endl;
					//n_color = colored_vertices.size() + 1;
					/*for (size_t i = 0; i < colored_vertices.size(); i++)
					{
						for (size_t j = 0; j < colored_vertices[i].size(); j++)
						{
							size_t local_id = colored_vertices[i][j];
							for (size_t k = 0; k < face_clusters[local_id].size(); k++)
							{
								face_color[face_clusters[local_id][k]] = i + 1;
							}
						}
					}*/
					for (size_t i = 0; i < face_clusters.size(); i++)
					{
						size_t cc = cluster_color[i];
						for (size_t j = 0; j < face_clusters[i].size(); j++)
						{
							face_color[face_clusters[i][j]] = cc + 1;
						}
					}
					
					update_tree_color(cluster_color, tree);

				}

				std::string tree_str = convert_tree_to_string<size_t>(tree);
				save_conf_file((output_prefix + "_csg.conf").c_str(), tree_str, flag_convex);
			}
			
			//sampling on triangles

			std::vector<double> tri_mean_curvature_normalize(tri_verts.size(), 0.0);
			if (result.count("cot"))
			{
				//compute tri_mean_curvature
				std::vector<double> vert_curvature;
				compute_vert_mean_curvature(vert_pos, tri_verts, vert_curvature);
				for (size_t i = 0; i < tri_verts.size(); i++)
				{
					for (size_t j = 0; j < 3; j++)
					{
						tri_mean_curvature_normalize[i] += vert_curvature[tri_verts[i][j]];
					}
				}
				//auto [minc, maxc] = std::minmax_element(begin(tri_mean_curvature_normalize), end(tri_mean_curvature_normalize));
				auto maxele = std::max_element(begin(tri_mean_curvature_normalize), end(tri_mean_curvature_normalize));
				auto minele = std::min_element(begin(tri_mean_curvature_normalize), end(tri_mean_curvature_normalize));
				std::cout << "min curvature: " << *minele << " max curvature: " << *maxele << std::endl;
				double diff = *maxele - *minele;
				for (size_t i = 0; i < tri_verts.size(); i++)
				{
					tri_mean_curvature_normalize[i] = (tri_mean_curvature_normalize[i] - *minele) / diff;
					//tri_mean_curvature_normalize[i] = (*maxele -  tri_mean_curvature_normalize[i]) / diff;
				}
			}

			for (size_t i = 0; i < tri_area.size(); i++)
			{
				tri_area[i] = std::abs(compute_tri_area<double>(vert_pos[tri_verts[i][0]], vert_pos[tri_verts[i][1]], vert_pos[tri_verts[i][2]]));
				if (result.count("cot"))
				{
					tri_area[i] = tri_area[i] * exp(5 * tri_mean_curvature_normalize[i]);
					//tri_area[i] = tri_area[i] * exp(5 * tri_mean_curvature_normalize[i]);
					
				}
				total_tri_area += tri_area[i];
			}
			std::vector<double> tri_bound(tri_area.size() + 1, 0.0);
			for (size_t i = 0; i < tri_area.size(); i++)
			{
				tri_bound[i + 1] = tri_bound[i] + tri_area[i] / total_tri_area;
			}

			//first select each point on each face, then sample the rest points randomly
			for (size_t i = 0; i < mesh.get_num_of_faces(); i++)
			{
				int fid = i;
				double s = unif_dist(e2);
				double t = unif_dist(e2);
				if (s + t > 1)
				{
					s = 1 - s;
					t = 1 - t;
				}
				TinyVector<double, 3> facenormal = mesh.get_faces_list()->at(fid)->normal;
				if (result.count("s"))
					sample_pts.push_back((1.0 - s - t) * vert_pos[tri_verts[fid][0]] + s * vert_pos[tri_verts[fid][1]] + t * vert_pos[tri_verts[fid][2]] + facenormal * normal_dist(e2));
				else
					sample_pts.push_back((1.0 - s - t) * vert_pos[tri_verts[fid][0]] + s * vert_pos[tri_verts[fid][1]] + t * vert_pos[tri_verts[fid][2]]);
				
				//sample mask to be added
				if (result.count("sn"))
				{
					sample_pt_normals.push_back(perturb_normal(facenormal, angle_unif_dist(e2), angle_unif_dist(e2)));
				}
				else
				{
					sample_pt_normals.push_back(facenormal);
				}
				assert(face_color[fid] <= n_color);
				sample_mask[n_feature_sample + i] = face_color[fid];
			}


			for (size_t i = mesh.get_num_of_faces(); i < n_nonfeature_sample; i++)
			{
				double u = unif_dist(e2);
				auto iter = std::upper_bound(tri_bound.begin(), tri_bound.end(), u);
				int fid = (int)std::distance(tri_bound.begin(), iter);
				assert(fid != tri_verts.size() + 1);
				fid = std::max(0, fid - 1);
				//sample
				//int id0 = ungrouped_features[fid].first;
				//int id1 = ungrouped_features[fid].second;
				double s = unif_dist(e2);
				double t = unif_dist(e2);
				if (s + t > 1)
				{
					s = 1 - s;
					t = 1 - t;
				}
				TinyVector<double, 3> facenormal = mesh.get_faces_list()->at(fid)->normal;
				if (result.count("s"))
					sample_pts.push_back((1.0 - s - t) * vert_pos[tri_verts[fid][0]] + s * vert_pos[tri_verts[fid][1]] + t * vert_pos[tri_verts[fid][2]] + facenormal * normal_dist(e2));
				else
					sample_pts.push_back((1.0 - s - t) * vert_pos[tri_verts[fid][0]] + s * vert_pos[tri_verts[fid][1]] + t * vert_pos[tri_verts[fid][2]]);
				//sample_pt_normals.push_back(facenormal);
				if (result.count("sn"))
				{
					sample_pt_normals.push_back(perturb_normal(facenormal, angle_unif_dist(e2), angle_unif_dist(e2)));
				}
				else
				{
					sample_pt_normals.push_back(facenormal);
				}
				//sample mask to be added
				assert(face_color[fid] <= n_color);
				sample_mask[n_feature_sample + i] = face_color[fid];

				//sample_pts.push_back(s * mesh.get_vertices_list()->at(id0)->pos + (1.0 - s) * mesh.get_vertices_list()->at(id1)->pos);
				//sample_pt_normals.push_back(TinyVector<double, 3>(1.0, 0.0, 0.0));
			}

			std::ofstream outputsamples(outputfile.c_str());
			
			for (size_t i = 0; i < sample_pts.size(); i++)
			{
				outputsamples << sample_pts[i] << " " << sample_pt_normals[i] << std::endl;
			}

			outputsamples.close();
			
			
			//output ply
			//double u = unif_dist(e2);
			std::vector<std::array<double, 3>> branch_color;
			for (size_t i = 0; i < n_color; i++)
			{
				branch_color.push_back(std::array<double, 3>{ {unif_dist(e2), unif_dist(e2), unif_dist(e2)}});
			}

			std::vector<std::array<double, 3>> meshVertexPositions;
			std::vector<std::array<double, 3>> meshVertexColors;

			for (size_t i = 0; i < sample_pts.size(); i++)
			{
				meshVertexPositions.push_back(std::array<double, 3>{ {sample_pts[i][0], sample_pts[i][1], sample_pts[i][2]}});
				meshVertexColors.push_back(branch_color[sample_mask[i]]);
			}


			// Create an empty object
			happly::PLYData plyOut;

			// Add mesh data (elements are created automatically)
			plyOut.addVertexPositions(meshVertexPositions);
			plyOut.addVertexColors(meshVertexColors);

			// Write the object to file
			plyOut.write(output_prefix + "_patch_"+ std::to_string(n_color - 1) +".ply", happly::DataFormat::ASCII);

			_mkdir(output_prefix.c_str());
			for (size_t i = 0; i < n_color; i++)
			{
				std::ofstream ofs((output_prefix + "\\" +  std::to_string(i) + ".xyz"));
				for (size_t j = 0; j < sample_pts.size(); j++)
				{
					if (sample_mask[j] == i)
					{
						ofs << sample_pts[j] << " " << sample_pt_normals[j] << std::endl;
					}
				}
				ofs.close();
			}

			std::ofstream outputmask(outputmaskfile.c_str());
			for (size_t i = 0; i < sample_mask.size(); i++)
			{
				outputmask << sample_mask[i] << std::endl;
			}

			outputmask.close();
			return 1;
		}

	}
	catch (const cxxopts::OptionException& e)
	{
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	}
	return 0;
	
	//int process_mode = std::stoi(argv[1]);
	//if (process_mode == 0 && argc != 4)
	//{
	//	std::cout << "Input format: " << argv[0] << " process_mode input.obj output.obj" << std::endl;
	//	return 1;
	//}

	//if (process_mode == 1 && argc != 7 )
	//{
	//	std::cout << "Input format: " << argv[0] << " process_mode input.obj input.fea output.points/xyz output_mask.txt [optional]n_nonfeature_sample" << std::endl;
	//}
	//


	//Mesh3d mesh;
	//mesh.load_obj(argv[2]);
	//std::cout << "verts: " << mesh.get_vertices_list()->size() << " face:  " << mesh.get_faces_list()->size() << std::endl;
	//std::cout << "process mode: " << process_mode << std::endl;
	//
	//if (process_mode == 0)
	//{
	//	std::cout << "enter normalizing" << std::endl;
	//	////normalization part begin
	//	//[-0.9, 9]^3
	//	std::vector<TinyVector<double, 3>> pts_nl(mesh.get_vertices_list()->size());
	//	double max_range = mesh.xmax - mesh.xmin;
	//	max_range = max_range < (mesh.ymax - mesh.ymin) ? (mesh.ymax - mesh.ymin) : max_range;
	//	max_range = max_range < (mesh.zmax - mesh.zmin) ? (mesh.zmax - mesh.zmin) : max_range;

	//	double xcenter = (mesh.xmin + mesh.xmax) / 2;
	//	double ycenter = (mesh.ymin + mesh.ymax) / 2;
	//	double zcenter = (mesh.zmin + mesh.zmax) / 2;
	//	std::cout << "center " << xcenter << " " << ycenter << " " << zcenter << std::endl;

	//	for (size_t i = 0; i < mesh.get_vertices_list()->size(); i++)
	//	{
	//		mesh.get_vertices_list()->at(i)->pos[0] = (mesh.get_vertices_list()->at(i)->pos[0] - xcenter) / max_range * 1.8;
	//		mesh.get_vertices_list()->at(i)->pos[1] = (mesh.get_vertices_list()->at(i)->pos[1] - ycenter) / max_range * 1.8;
	//		mesh.get_vertices_list()->at(i)->pos[2] = (mesh.get_vertices_list()->at(i)->pos[2] - zcenter) / max_range * 1.8;
	//	}
	//	
	//	//output mesh
	//	mesh.write_obj(argv[2]);
	//	return 1;
	//	////normalization part end
	//}
	//

	////feature transfer part
	////currently all normals are set to 1.0 0.0 0.0
	//std::ifstream inputfeature(argv[2]);
	//std::ofstream outputfile(argv[3]);
	//std::ofstream outputmask(argv[4]);
	//int non_feature_samples = atoi(argv[5]);

	////std::vector<double> coordx, coordy, coordz;
	//std::vector<int> faces;
	///*std::string str;
	//inputfile >> str;
	//inputfile >> str;*/

	//std::vector<std::pair<int, int>> ungrouped_features;
	//std::vector<std::vector<std::pair<int, int>>> grouped_features;
	//load_feature_file(argv[2], ungrouped_features, grouped_features);

	////sample 10 pts on each feature segments
	//int n_seg = 10;
	//for (size_t i = 0; i < ungrouped_features.size(); i++)
	//{
	//	int id0 = ungrouped_features[i].first;
	//	int id1 = ungrouped_features[i].second;
	//	double sx = (mesh.get_vertices_list()->at(id1)->pos[0] - mesh.get_vertices_list()->at(id0)->pos[0]) / n_seg;
	//	double sy = (mesh.get_vertices_list()->at(id1)->pos[1] - mesh.get_vertices_list()->at(id0)->pos[1]) / n_seg;
	//	double sz = (mesh.get_vertices_list()->at(id1)->pos[2] - mesh.get_vertices_list()->at(id0)->pos[2]) / n_seg;
	//	//point normal as average normal of neighboring faces
	//	HE_edge<double>* te = NULL;
	//	HE_edge<double>* ve_begin = mesh.get_vertices_list()->at(id0)->edge;
	//	assert(ve_begin->pair->vert->id == id0);
	//	HE_edge<double>* ve_iter = ve_begin;
	//	do
	//	{
	//		int next_id = ve_iter->vert->id;
	//		if (next_id == id1)
	//		{
	//			te = ve_iter;
	//			break;
	//		}
	//		ve_iter = ve_iter->pair->next;
	//	} while (ve_iter != ve_begin);

	//	//std::cout << "edge id: " << te->id << std::endl;
	//	assert(te != NULL);
	//	TinyVector<double, 3> avg_normal = te->face->normal + te->pair->face->normal;
	//	avg_normal.Normalize();
	//	for (size_t j = 0; j < n_seg; j++)
	//	{
	//		outputfile << mesh.get_vertices_list()->at(id0)->pos[0] + sx * j << " " << mesh.get_vertices_list()->at(id0)->pos[1] + sy * j << " " << mesh.get_vertices_list()->at(id0)->pos[2] + sz * j << " " << avg_normal << std::endl;
	//	}
	//}

	//outputfile.close();
	//
	//for (size_t i = 0; i < non_feature_samples; i++)
	//{
	//	outputmask << "0\n";
	//}
	//for (size_t i = 0; i < ungrouped_features.size() * n_seg; i++)
	//{
	//	outputmask << "1\n";
	//}
	//outputmask.close();

	//system("pause");

	return 0;
}