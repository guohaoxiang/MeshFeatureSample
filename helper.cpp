#include "helper.h"
#include <algorithm>
#include <queue>
#include <igl/gaussian_curvature.h>
#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/readOFF.h>

void greedy_graph_coloring(const size_t num_vertices, const std::vector<std::set<size_t>>& edges, std::vector<std::vector<size_t>>& colored_vertices)
{
	std::vector<int> result(num_vertices, -1);
	result[0] = 0;
	std::vector<bool> available(num_vertices, false);

	int max_color = 0;
	for (size_t u = 1; u < num_vertices; u++)
	{
		for (typename std::set<size_t>::iterator iter = edges[u].begin(); iter != edges[u].end(); iter++)
		{
			if (result[*iter] != -1)
				available[result[*iter]] = true;
		}
		int cr;
		for (cr = 0; cr < (int)num_vertices; cr++)
		{
			if (available[cr] == false)
				break;
		}

		result[u] = cr;
		max_color = std::max(cr, max_color);

		for (typename std::set<size_t>::iterator iter = edges[u].begin(); iter != edges[u].end(); iter++)
		{
			if (result[*iter] != -1)
				available[result[*iter]] = false;
		}
	}
	colored_vertices.resize(max_color + 1);
	for (size_t i = 0; i < num_vertices; i++)
	{
		colored_vertices[result[i]].push_back(i);
	}

}

void get_graph_component(const std::vector<std::set<size_t>>& edges, std::vector<std::set<size_t>>& components)
{
	//edges id starting from 0
	components.clear();
	std::vector<int> vert2color(edges.size(), -1);
	int cid = 0;
	int startid = -1;
	for (size_t i = 0; i < edges.size(); i++)
	{
		if (!edges[i].empty())
		{
			startid = i;
			break;
		}
	}
	assert(startid != -1);
	while (startid != -1)
	{
		std::set<size_t> onecluser;
		std::queue<size_t> q;
		q.push(startid);
		vert2color[startid] = cid;
		while (!q.empty())
		{
			size_t front = q.front();
			q.pop();
			onecluser.insert(front);
			for (auto v : edges[front])
			{
				if (vert2color[v] == -1)
				{
					q.push(v);
					vert2color[v] = cid;
				}
			}
		}
		components.push_back(onecluser);
		cid++;
		startid = -1;
		for (size_t i = 0; i < edges.size(); i++)
		{
			if (!edges[i].empty() && vert2color[i] == -1)
			{
				startid = i;
				break;
			}
		}
	}
}

void compute_vert_mean_curvature(const std::vector<TinyVector<double, 3>>& pos, const std::vector<TinyVector<size_t, 3>>& faces, std::vector<double>& curvature)
{
	using namespace Eigen;
	using namespace std;
	MatrixXd V;
	MatrixXi F;
	curvature.clear();
	curvature.resize(pos.size(), 0.0);
	//igl::readOFF(TUTORIAL_SHARED_PATH "/bumpy.off",V,F);
	//igl::readOFF("E:\\code\\PaperCode\\VolumeMeshProcessing\\x64\\Release\\AllModels\\fandisk\\fandisk.off", V, F);
	
	
	//remove dup verts
	std::vector<size_t> o2n(pos.size(), (size_t)-1), n2o;
	size_t count = 0;
	for (size_t i = 0; i < faces.size(); i++)
	{
		for (size_t j = 0; j < 3; j++)
		{
			size_t id = faces[i][j];
			if (o2n[id] == (size_t)-1)
			{
				o2n[id] = count;
				count = count + 1;
				n2o.push_back(id);
			}
		}
	}

	std::vector<std::vector<double>> vV(n2o.size(), std::vector<double>(3));
	std::vector<std::vector<int>> vF(faces.size(), std::vector<int>(3));
	for (size_t i = 0; i < n2o.size(); i++)
	{
		//V(3 *i, 1) = pos[i][0];
		/*V(3 * i + 1) = pos[i][1];
		V(3 * i + 2) = pos[i][2];*/
		for (size_t j = 0; j < 3; j++)
		{
			vV[i][j] = pos[n2o[i]][j];
		}
	}
	for (size_t i = 0; i < faces.size(); i++)
	{
		/*F(3 * i) = faces[i][0];
		F(3 * i + 1) = faces[i][1];
		F(3 * i + 2) = faces[i][2];*/
		for (size_t j = 0; j < 3; j++)
		{
			vF[i][j] = o2n[faces[i][j]];
		}
	}
	igl::list_to_matrix(vV, V);
	igl::list_to_matrix(vF, F);
	std::cout << "V size: " << V.size() << std::endl;
	std::cout << "F size: " << F.size() << std::endl;
	std::cout << "pos & faces size: " << pos.size() << " " << faces.size() << std::endl;


	SparseMatrix<double> M, Minv;
	MatrixXd HN;
	SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
	igl::invert_diag(M, Minv);
	HN = -Minv * (L * V);
	VectorXd H = HN.rowwise().norm(); //up to sign
	std::cout << "H min max: " << H.minCoeff() << " " << H.maxCoeff() << std::endl;
	for (size_t i = 0; i < n2o.size(); i++)
	{
		curvature[n2o[i]] = H[i];
	}
}

void update_tree_color(const std::vector<size_t>& new_color, TreeNode<size_t>* t)
{
	std::set<size_t> new_keys;
	for (auto v : t->keys)
	{
		new_keys.insert(new_color[v]);
	}
	t->keys = new_keys;
	for (size_t i = 0; i < t->children.size(); i++)
	{
		update_tree_color(new_color, t->children[i]);
	}

}



bool get_tree_from_convex_graph(const std::vector<std::set<size_t>> &graph , const std::map<std::pair<size_t, size_t>, int>& flag_fpconvex, bool flag_convex_bool, TreeNode<size_t>* tn, int layer)
{
	//set tn
	int target_convex = 2 - (int)flag_convex_bool;
	std::set<size_t> counter_nodes;
	for (size_t i = 0; i < graph.size(); i++)
	{
		for (auto v : graph[i])
		{
			if (v > i)
			{
				auto it = flag_fpconvex.find(std::pair<size_t, size_t>(i, v));
				assert(it != flag_fpconvex.end());
				if (it->second != 0 && it->second != target_convex)
				{
					counter_nodes.insert(i);
					counter_nodes.insert(v);
				}
			}
			
		}
	}

	if (!counter_nodes.empty())
	{
		std::vector<std::set<size_t>> counter_graph(graph.size()); //graph made of counter nodes
		for (auto cv : counter_nodes)
		{
			for (auto cvn : graph[cv])
			{
				if (std::find(counter_nodes.begin(), counter_nodes.end(), cvn) != counter_nodes.end())
				{
					auto it = flag_fpconvex.find(std::pair<size_t, size_t>(std::min(cv, cvn), std::max(cv, cvn)));
					assert(it != flag_fpconvex.end());
					if (it->second != target_convex) //including smooth one
						counter_graph[cv].insert(cvn);
				}
			}
		}
		std::vector<std::set<size_t>> counter_clusters;
		get_graph_component(counter_graph, counter_clusters);
		//set tn.children
		for (size_t i = 0; i < counter_clusters.size(); i++)
		{
			std::vector<std::set<size_t>> subgraph(graph.size());
			for (auto v : counter_clusters[i])
			{
				for (auto vn : graph[v])
				{
					if (std::find(counter_clusters[i].begin(), counter_clusters[i].end(), vn) != counter_clusters[i].end())
					{
						subgraph[v].insert(vn);
					}
				}
			}
			TreeNode<size_t>* child = new TreeNode<size_t>;
			if (layer == 10)
			{
				std::cout << "Layers over 10!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
				//exit(EXIT_FAILURE);
				return false;
			}
			bool tmp_flag = get_tree_from_convex_graph(subgraph, flag_fpconvex, !flag_convex_bool, child, layer + 1);
			if (!tmp_flag)
			{
				return false;
			}
			tn->children.push_back(child);
		}
	}

	//set nodes that are not counter nodes to tn.keys
	for (size_t i = 0; i < graph.size(); i++)
	{
		if (!graph[i].empty() && std::find(counter_nodes.begin(), counter_nodes.end(), i) == counter_nodes.end())
		{
			tn->keys.insert(i);
		}
	}
	return true;
}

using namespace MeshLib;
void sort_grouped_features(Mesh3d* m, std::vector<std::vector<int>>& grouped_features)
{
	std::vector<std::vector<int>> grouped_features_new;
	for (size_t i = 0; i < grouped_features.size(); i++)
	{
		std::map<int, std::vector<int>> v2hes;
		std::map<int, std::vector<int>> v2vs;
		for (size_t j = 0; j < grouped_features[i].size(); j++)
		{
			int heid = grouped_features[i][j];
			int vid[2];
			vid[0] = m->get_edges_list()->at(heid)->vert->id;
			vid[1] = m->get_edges_list()->at(heid)->pair->vert->id;
			for (size_t k = 0; k < 2; k++)
			{
				if (v2hes.find(vid[k]) == v2hes.end())
				{
					std::vector<int> hes, vs;
					hes.push_back(heid);
					vs.push_back(vid[1 - k]);
					/*v2hes[vid[k]] = std::vector<int>(heid);
					v2vs[vid[k]] = std::vector<int>(vid[1 - k]);*/
					v2hes[vid[k]] = hes;
					v2vs[vid[k]] = vs;
				}
				else
				{
					v2hes[vid[k]].push_back(heid);
					v2vs[vid[k]].push_back(vid[1 - k]);
				}
			}
		}

		int cur_edge = -1;
		int cur_vert = -1;
		int prev_edge = -1;
		for (auto it : v2hes)
		{
			if (it.second.size() == 1)
			{
				cur_vert = it.first;
			}
		}
		//assert(cur_edge != -1);
		if (cur_vert == -1)
		{
			//circular case
			cur_vert = v2hes.begin()->first;
		}
		std::vector<int> one_group;
		//from start to end
		std::queue<int> q;
		std::vector<bool> edge_color(m->get_num_of_edges(), false);
		while (true)
		{
			//get cur_edge
			cur_edge = -1;
			for (auto e : v2hes[cur_vert])
			{
				if (e != prev_edge && edge_color[e] == false)
				{
					cur_edge = e;
					break;
				}
			}
			if (cur_edge == -1)
				break;

			edge_color[cur_edge] = true;
			//get next_vert
			assert(cur_vert == m->get_edges_list()->at(cur_edge)->vert->id || cur_vert == m->get_edges_list()->at(cur_edge)->pair->vert->id);

			int next_vert = m->get_edges_list()->at(cur_edge)->vert->id + m->get_edges_list()->at(cur_edge)->pair->vert->id - cur_vert;
			
			if (next_vert == m->get_edges_list()->at(cur_edge)->vert->id)
			{
				one_group.push_back(cur_edge);
			}
			else
			{
				one_group.push_back(m->get_edges_list()->at(cur_edge)->pair->id);
			}

			cur_vert = next_vert;
			prev_edge = cur_edge;
			
			
		}
		
		assert(one_group.size() == grouped_features[i].size());
		grouped_features_new.push_back(one_group);
	}
	grouped_features = grouped_features_new;
}

int cluster_mesh_faces(Mesh3d* m, const std::vector<bool>& he_feature_flag, std::vector<std::vector<int>> &grouped_features, int cluster_begin_id, std::vector<int>& face2cluster, std::vector<std::pair<int, int>>& feature_twoface_colors)
{
	//return n_cluster + cluster_start_id
	face2cluster.clear();
	feature_twoface_colors.clear();
	face2cluster.resize(m->get_num_of_faces(), -1);
	std::vector<TinyVector<double, 3>> face_centers(m->get_num_of_faces(), TinyVector<double, 3>(0,0,0));
	for (size_t i = 0; i < m->get_num_of_faces(); i++)
	{
		HE_edge<double>* begin_edge = m->get_faces_list()->at(i)->edge;
		HE_edge<double>* edge = begin_edge;
		do
		{
			face_centers[i] = face_centers[i] + edge->vert->pos;
			edge = edge->next;
		} while (edge != begin_edge);
		face_centers[i] = face_centers[i] / 3.0;
	}


	class face_cluster
	{
	public:
		face_cluster(int faceid = 0,  TinyVector<double, 3> ori = TinyVector<double, 3>(0,0,0), TinyVector<double, 3> cur = TinyVector<double, 3>(0, 0, 0))
			:fid(faceid), ori_face_center(ori)
		{
			e = (cur - ori).Length();
			//e = energy;
		}
	public:
		int fid;
		TinyVector<double, 3> ori_face_center;
		double e;
	};
	
	class compare_face_cluster
	{
	public:
		compare_face_cluster()
		{
		}
		bool operator()(const face_cluster& fc1, const face_cluster& fc2) const
		{
			return fc1.e > fc2.e;
		}
	};
	
	std::vector<double> energy(m->get_num_of_faces(), DBL_MAX);
	std::priority_queue<face_cluster, std::deque<face_cluster>, compare_face_cluster> m_queue;
	sort_grouped_features(m, grouped_features);
	std::vector<bool> flag_face_unchangable(m->get_num_of_faces(), false);
	for (size_t i = 0; i < grouped_features.size(); i++)
	{
		int color1 = cluster_begin_id;
		int color2 = cluster_begin_id + 1;
		feature_twoface_colors.push_back(std::pair<int, int>(color1, color2));
		for (size_t j = 0; j < grouped_features[i].size(); j++)
		{
			int heid = grouped_features[i][j];
			int fid1 = m->get_edges_list()->at(heid)->face->id;
			int fid2 = m->get_edges_list()->at(heid)->pair->face->id;

			if (face2cluster[fid1] == -1)
			{
				face2cluster[fid1] = color1;
				m_queue.push(face_cluster(fid1, face_centers[fid1], face_centers[fid1]));
				energy[fid1] = 0;
				/*energy[fid1] = 1 - std::abs(m->get_faces_list()->at(fid1)->normal.Dot(m->get_faces_list()->at(fid2)->normal));
				m_queue.push(face_cluster(fid1, face_centers[fid1], face_centers[fid1], energy[fid1]));*/
				flag_face_unchangable[fid1] = true;
			}
			if (face2cluster[fid2] == -1)
			{
				face2cluster[fid2] = color2;
				m_queue.push(face_cluster(fid2, face_centers[fid2], face_centers[fid2]));
				energy[fid2] = 0;
				//energy[fid2] = energy[fid1];
				//m_queue.push(face_cluster(fid2, face_centers[fid2], face_centers[fid2], energy[fid2]));

				flag_face_unchangable[fid2] = true;
			}
		}
		cluster_begin_id = cluster_begin_id + 2;
	}

	while (!m_queue.empty())
	{
		face_cluster fc = m_queue.top();
		m_queue.pop();
		HE_edge<double>* begin_edge = m->get_faces_list()->at(fc.fid)->edge;
		HE_edge<double>* edge = begin_edge;
		do
		{
			if (he_feature_flag[edge->id])
			{
				edge = edge->next;
				continue;
			}
			assert(fc.fid == edge->face->id || fc.fid == edge->pair->face->id);
			int otherface = edge->face->id + edge->pair->face->id - fc.fid;
			if (!flag_face_unchangable[otherface])
			{
				double tmp_energy = (face_centers[otherface] - fc.ori_face_center).Length();
				//double tmp_energy = 1.0 - std::abs(edge->face->normal.Dot(edge->pair->face->normal));
				if (face2cluster[otherface] == -1)
				{
					face2cluster[otherface] = face2cluster[fc.fid];
					m_queue.push(face_cluster(otherface, fc.ori_face_center, face_centers[otherface]));
					energy[otherface] = tmp_energy;
					//m_queue.push(face_cluster(otherface, fc.ori_face_center, face_centers[otherface], tmp_energy));
				}
				else
				{
					//colored
					if (tmp_energy < energy[otherface])
					{
						face2cluster[otherface] = face2cluster[fc.fid];
						//m_queue.push(face_cluster(otherface, tmp_energy));
						m_queue.push(face_cluster(otherface, fc.ori_face_center, face_centers[otherface]));
						energy[otherface] = tmp_energy;
						//m_queue.push(face_cluster(otherface, fc.ori_face_center, face_centers[otherface], tmp_energy));

					}
				}
			}
			edge = edge->next;
		} while (edge != begin_edge);
	}
	return cluster_begin_id;
}

void get_cluster_from_coloring(const std::vector<int>& face_color, int color_start_id, std::vector<std::vector<int>>& face_clusters)
{
	int max_color = -1;
	for (auto c : face_color)
	{
		assert(c >= color_start_id);
		if (max_color < c)
		{
			max_color = c;
		}
	}
	face_clusters.clear();
	for (size_t i = color_start_id; i <= max_color; i++)
	{
		std::vector<int> one_cluster;
		for (size_t j = 0; j < face_color.size(); j++)
		{
			if (face_color[j] == i)
			{
				one_cluster.push_back(j);
			}
		}
		face_clusters.push_back(one_cluster);
	}
}

int merge_clusters(Mesh3d* m, const std::vector<bool>& he_feature_flag, int cluster_begin_id, int n_cluster, const std::vector<std::pair<int, int>>& feature_twoface_colors, std::vector<int>& face_color)
{
	//return num_clusters + cluster_begin_id
	std::vector<int> feature_pair_mapping(n_cluster, -1);
	for (auto& p : feature_twoface_colors)
	{
		feature_pair_mapping[p.first - 1] = p.second - 1;
		feature_pair_mapping[p.second - 1] = p.first - 1;
	}
	std::vector<int> cluster_o2n(n_cluster, -1); //starting from 0->0
	std::vector<std::set<int>> con(n_cluster); //split by features
	int cluster_id = 0;

	for (size_t i = 0; i < m->get_num_of_edges(); i++)
	{
		if (he_feature_flag[i])
			continue;
		int color1 = face_color[m->get_edges_list()->at(i)->face->id] - 1;
		int color2 = face_color[m->get_edges_list()->at(i)->pair->face->id] - 1;
		if (color1 != color2)
		{
			con[color1].insert(color2);
			con[color2].insert(color1);
		}
	}
	
	for (size_t i = 0; i < n_cluster; i++)
	{
		if (cluster_o2n[i] != -1)
			continue;
		std::queue<int> q;
		q.push(i);
		cluster_o2n[i] = cluster_id;
		std::set<int> mset;
		while (!q.empty())
		{
			int front = q.front();
			q.pop();
			for (auto nn : con[front])
			{
				if (cluster_o2n[nn] != -1)
					continue;
				int pair_cluster = feature_pair_mapping[nn];
				if (cluster_o2n[pair_cluster] == cluster_id)
					continue;
				cluster_o2n[nn] = cluster_id;
				q.push(nn);
			}
		}
		cluster_id = cluster_id + 1;
	}

	for (size_t i = 0; i < face_color.size(); i++)
	{
		face_color[i] = cluster_o2n[face_color[i] - 1] + 1;
	}

	return cluster_id + 1;
}

void get_grouped_edges(Mesh3d& mesh, const std::vector<bool>& he_feature_flag, const std::vector<std::vector<int>>& feature_v2he, std::vector<std::vector<int>>& grouped_features, std::vector<int>& he2gid)
{
	grouped_features.clear();
	he2gid.clear();
	he2gid.resize(mesh.get_num_of_edges(), -1);
	int groupid = 0;
	while (true)
	{
		int first_he = -1;
		for (size_t i = 0; i < he_feature_flag.size(); i++)
		{
			if (he_feature_flag[i] && he2gid[i] == -1)
			{
				first_he = i;
				break;
			}
		}
		if (first_he == -1)
		{
			break;
		}

		std::vector<int> onegroup;
		std::queue<size_t> q;
		q.push(first_he);
		he2gid[first_he] = groupid;
		he2gid[mesh.get_edges_list()->at(first_he)->pair->id] = groupid;
		while (!q.empty())
		{
			size_t curhe = q.front();
			onegroup.push_back(curhe);
			q.pop();
			size_t curhe_pair = mesh.get_edges_list()->at(curhe)->pair->id;
			std::vector<size_t> twoverts;
			twoverts.push_back(mesh.get_edges_list()->at(curhe)->vert->id);
			twoverts.push_back(mesh.get_edges_list()->at(curhe)->pair->vert->id);
			std::vector<size_t> twoedges;
			twoedges.push_back(curhe_pair);
			twoedges.push_back(curhe);
			for (size_t i = 0; i < 2; i++)
			{
				//only accept points with degree 2
				if (feature_v2he[twoverts[i]].size() == 2)
				{
					assert(twoedges[i] == feature_v2he[twoverts[i]][0] || twoedges[i] == feature_v2he[twoverts[i]][1]);
					size_t other_he = feature_v2he[twoverts[i]][0] + feature_v2he[twoverts[i]][1] - twoedges[i];
					if (he2gid[other_he] == -1)
					{
						q.push(other_he);
						he2gid[other_he] = groupid;
						he2gid[mesh.get_edges_list()->at(other_he)->pair->id] = groupid;
					}
				}
			}
		}
		grouped_features.push_back(onegroup);
		groupid++;
	}
}