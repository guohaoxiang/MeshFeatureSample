#include <iostream>
#include <fstream>
#include "Mesh3D.h"
#include "TinyVector.h"
#include "happly/happly.h"
#include "cxxopts.hpp"

#define SHORTEST_EDGE_LENGTH 1e-6


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

void select_features_by_angle(Mesh3d* m, double th_angle, std::vector<std::pair<int, int>> &ungrouped_features, std::vector<int> &ungrouped_feature_eid)
{
	ungrouped_features.clear();
	ungrouped_feature_eid.clear();
	const double angle = th_angle * 3.1415926535897932384626433832795 / 180;
	std::vector<bool> he_flag(m->get_num_of_edges(), false);
	for (size_t i = 0; i < m->get_num_of_edges(); i++)
	{
		HE_edge<double>* he = m->get_edge(i);
		if (!he_flag[he->id])
		{
			if (!m->is_on_boundary(he) && acos(he->face->normal.Dot(he->pair->face->normal)) > angle)
			{
				he_flag[he->id] = true;
				he_flag[he->pair->id] = true;
				ungrouped_features.push_back(std::pair<int, int>(he->vert->id, he->pair->vert->id));
				ungrouped_feature_eid.push_back(he->id);
			}
		}
	}
}

void output_pts_xyz(const std::vector<TinyVector<double, 3>>& pts, const char* filename)
{
	std::ofstream ofs(filename);
	for (size_t i = 0; i < pts.size(); i++)
	{
		ofs << pts[i] << std::endl;
	}

	ofs.close();
}

int main(int argc, char** argv)
{
	try
	{
		cxxopts::Options options("SimpleSample", "Sampling feature points from points on mesh (author: Haoxiang Guo, Email: guohaoxiangxiang@gmail.com)");
		options
			.positional_help("[optional args]")
			.show_positional_help()
			.allow_unrecognised_options()
			.add_options()
			("i,input", "input mesh (obj/off/ply format)", cxxopts::value<std::string>())
			("f,feature", "input feature file (fea format)", cxxopts::value<std::string>())
			("o,output", "output points (ptangle format)", cxxopts::value<std::string>())
			("a", "angle threshold in degree for detecting features, default(30)", cxxopts::value<double>())
			("s", "length of line segment for sampling, default(4e-3)", cxxopts::value<double>())
			("h,help", "print help");

		auto result = options.parse(argc, argv);
		if (result.count("help"))
		{
			std::cout << options.help({ "", "Group" }) << std::endl;
			exit(0);
		}
		bool flag_feature_flag = true;
		double len_seg = 4e-3;
		double th_angle = 30;
		if (result.count("s"))
		{
			len_seg = result["s"].as<double>();
		}
		if (result.count("a"))
		{
			//angle is used only when no feature files are given
			th_angle = result["a"].as<double>();
		}

		auto& inputfile = result["i"].as<std::string>();
		auto& outputfile = result["o"].as<std::string>();

		int last_dot = (int)outputfile.find_last_of(".");
		auto output_prefix = outputfile.substr(0, last_dot);
		
		std::string inputext = GetFileExtension(inputfile);
		Mesh3d mesh;
		if (inputext == "obj")
			mesh.load_obj(inputfile.c_str());
		else if (inputext == "off")
			mesh.load_off(inputfile.c_str());
		else if (inputext == "ply")
		{
			happly::PLYData plyIn(inputfile);
			std::vector<std::array<double, 3>> vPos = plyIn.getVertexPositions();
			std::vector<std::vector<size_t>> fInd = plyIn.getFaceIndices<size_t>();
			mesh.load_mesh(vPos, fInd);

			//test code below
			//mesh.write_obj("test.obj");
		}
		std::cout << "verts: " << mesh.get_vertices_list()->size() << " face:  " << mesh.get_faces_list()->size() << std::endl;

		std::vector<std::pair<int, int>> ungrouped_features;
		std::vector<int> ungrouped_feature_eid;
		if (result.count("f"))
		{
			auto& inputfeaturefile = result["f"].as<std::string>();
			load_feature_file(inputfeaturefile.c_str(), ungrouped_features);
		}
		else
		{
			//select feature by angle
			select_features_by_angle(&mesh, th_angle, ungrouped_features, ungrouped_feature_eid);
		}
		
		//sample points
		std::vector<TinyVector<double, 3>> sample_pts;
		std::vector<double> pts_angle;
		//for (auto& pair : ungrouped_features)
		for (size_t i = 0; i < ungrouped_features.size(); i++)
		{
			auto pair = ungrouped_features[i];
			TinyVector<double, 3> pos0 = mesh.get_vertex(pair.first)->pos;
			TinyVector<double, 3> pos1 = mesh.get_vertex(pair.second)->pos;
			double len = (pos0 - pos1).Length();
			if (len > SHORTEST_EDGE_LENGTH)
			{
				//sample at least one points
				double angle = 0.0;
				if (ungrouped_feature_eid.empty())
				{
					bool flag = false;

					HE_edge<double>* begin_edge = mesh.get_vertex(pair.first)->edge;
					HE_edge<double>* edge = begin_edge;
					do
					{
						if (edge->vert->id == pair.second)
						{
							flag = true;
							angle = acos(edge->face->normal.Dot(edge->pair->face->normal)) * 180.0 / 3.1415926535897932384626433832795;
							break;
						}

						edge = edge->pair->next;
					} while (edge != begin_edge);

					assert(flag == true);
				}
				else
				{
					HE_edge<double>* edge = mesh.get_edge(ungrouped_feature_eid[i]);
					angle = acos(edge->face->normal.Dot(edge->pair->face->normal)) * 180.0 / 3.1415926535897932384626433832795;
				}
				
				if (len < len_seg / 2.0)
				{
					sample_pts.push_back((pos0 + pos1) / 2.0);
					pts_angle.push_back(angle);
				}
				else
				{
					int n_split = std::ceil((len - len_seg / 2.0) / len_seg);
					TinyVector<double, 3> vec =  (pos1 - pos0)/(pos1 - pos0).Length();
					for (size_t i = 0; i < n_split; i++)
					{
						double tmp_len = len_seg / 2.0 + i * len_seg;
						sample_pts.push_back(pos0 + tmp_len * vec);
						pts_angle.push_back(angle);
					}
				}
			}
		}
		assert(pts_angle.size() == sample_pts.size());

		//output_xyz
		output_pts_xyz(sample_pts, (output_prefix + ".xyz").c_str());
		
		std::ofstream ofs(outputfile);
		for (size_t i = 0; i < pts_angle.size(); i++)
		{
			ofs << sample_pts[i] << " " << pts_angle[i] << std::endl;
		}

		ofs.close();

		//smoothness term
	}
	catch (const cxxopts::OptionException& e)
	{
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	}

	return 0;
}