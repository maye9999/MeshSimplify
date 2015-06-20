#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <set>
#include <cassert>
#include <map>
#include <ctime>
#include "Matrix.h"
#include "Vector.h"

using namespace std;

const double threshold = 100;

struct Vertex;
struct Edge;
struct Face;
struct Node;

vector<Vertex> vertices_pool;
vector<Edge> edges_pool;
vector<Face> faces_pool;
set<Node> BST;


struct Vertex
{
	vector<int> faces_id;
	vector<int> edges_id;
	Vector pos;
	Matrix q;
	bool enable;
	//Vertex() : enable(true) {}
	Vertex(Vector p) : enable(true), pos(p) {}

	void updateQMatrix();
};

struct Edge
{
	int vertices_id[2];
	Vector v_mid;
	double err;
	bool enable;
	//Edge() : enable(true) {}
	Edge(int a, int b) : enable(true) {vertices_id[0] = min(a, b); vertices_id[1] = max(a, b);}

	bool operator==(const Edge& e) const
	{
		return ((vertices_id[0] == e.vertices_id[0]) && (vertices_id[1] == e.vertices_id[1]));
	}
	bool operator<(const Edge& e) const
	{
		if (vertices_id[0] != e.vertices_id[0])
			return vertices_id[0] < e.vertices_id[0];
		return vertices_id[1] < e.vertices_id[1];
	}
	void confine()
	{
		int a = vertices_id[0], b = vertices_id[1];
		vertices_id[0] = min(a, b);
		vertices_id[1] = max(a, b);
	}
	void updateNewPoint()
	{
		assert(enable);
		// For Present, just get the mid-point
		assert(vertices_pool[vertices_id[0]].enable);
		assert(vertices_pool[vertices_id[1]].enable);
		v_mid = 0.5 * (vertices_pool[vertices_id[0]].pos + vertices_pool[vertices_id[1]].pos);
	}

	// invoke updateNewPoint() First!
	void updateErr();

	bool checkValid()
	{
		return (vertices_id[0] != vertices_id[1]);
	}


};

struct Face
{
	int vertices_id[3];
	Vector norm;
	Matrix kp;
	bool enable;
	//Face() : enable(true) {}
	Face(int a, int b, int c) : enable(true) 
	{
		vertices_id[0] = min(a, min(b, c));
		vertices_id[2] = max(a, max(b, c));
		vertices_id[1] = a + b + c - vertices_id[0] - vertices_id[2];
		updateNormal();
	}
	bool operator==(const Face& f)
	{
		return vertices_id[0] == f.vertices_id[0] 
			&& vertices_id[1] == f.vertices_id[1] 
			&& vertices_id[2] == f.vertices_id[3];
	}
	void confine()
	{
		int a = vertices_id[0], b = vertices_id[1], c = vertices_id[2];
		vertices_id[0] = min(a, min(b, c));
		vertices_id[2] = max(a, max(b, c));
		vertices_id[1] = a + b + c - vertices_id[0] - vertices_id[2];
	}
	void updateKp()
	{
		assert(enable);
		kp.clear();
		Matrix::mul_fast(norm, norm, kp);
	}
	void updateNormal()
	{
		Vector va = vertices_pool[vertices_id[0]].pos;
		Vector vb = vertices_pool[vertices_id[1]].pos;
		Vector vc = vertices_pool[vertices_id[2]].pos;
		norm = Vector::normalize(Vector::cross(va - vb, vb - vc));
		norm[3] = -1.0 * (norm[0] * va[0] + norm[1] * va[1] + norm[2] * va[2]);
	}
	bool checkValid()
	{
		return (vertices_id[0] != vertices_id[1]
			&& vertices_id[0] != vertices_id[2]
			&& vertices_id[1] != vertices_id[2]);
	}


};

struct Node
{
	int edge_index;
	Node(int i) : edge_index(i) {}
	bool operator<(const Node& n) const
	{
		double e1 = edges_pool[edge_index].err;
		double e2 = edges_pool[n.edge_index].err;
		if(e1 != e2)
			return e1 < e2;
		return edge_index < n.edge_index;
	}
};


inline void Vertex::updateQMatrix()
{
	assert(enable);
	q.clear();
	for(int j = 0; j < faces_id.size(); ++j)
	{
		Face& f = faces_pool[faces_id[j]];
		if(!f.enable)
			continue;
		q += faces_pool[faces_id[j]].kp;
	}	
}

void Edge::updateErr()
{
	assert(enable);
	err = 0;
	for(int j = 0; j < 2; ++j)
		for(int i : vertices_pool[vertices_id[j]].faces_id)
			err += Matrix::mul(v_mid, faces_pool[i].kp) * v_mid;
}

bool loadOBJFile(const char* f_name)
{
	FILE* fp = fopen(f_name, "r");
	if(fp==NULL)
	{
		printf("Error: Loading %s failed.\n",f_name);
		return false;
	}
	char buf[256];
	int lineNumber = 0;
	while(fscanf(fp, "%s", buf) != EOF)
	{
		lineNumber ++;
		switch(buf[0])
		{
		case '#':				/* comment */
			/* eat up rest of line */
			fgets(buf, sizeof(buf), fp);
			break;
		case 'v':				/* v*/
			double x, y, z;
			if(fscanf(fp, "%lf %lf %lf", &x, &y, &z)==3)
			{
				Vector v(x, y, z);
				vertices_pool.push_back(Vertex(v));
			}
			else 
			{
				printf("Error: Wrong Number of Values(Should be 3). at Line %d\n",lineNumber);
				return false;
			}	
			break;
		case 'f':				/* face */
			{
				int a, b, c;
				if( fscanf(fp, "%d", &a) ==1 && 
					fscanf(fp, "%d", &b)  ==1 &&
					fscanf(fp, "%d", &c)  ==1 )
				{
					a--; b--; c--;
					Face f(a, b, c);
					faces_pool.push_back(f);
					vertices_pool[a].faces_id.push_back(faces_pool.size()-1);
					vertices_pool[b].faces_id.push_back(faces_pool.size()-1);
					vertices_pool[c].faces_id.push_back(faces_pool.size()-1);
					Edge e(a, b);
					edges_pool.push_back(e);
					Edge ee(b, c);
					edges_pool.push_back(ee);
					Edge eee(a, c);
					edges_pool.push_back(eee);
				}
				else
				{
					printf("Error: Wrong Face at Line %d\n",lineNumber);
					return false;
				}
			}
			break;
		default:
			/* eat up rest of line */
			fgets(buf, sizeof(buf), fp);
			break;
		}
	}
	sort(edges_pool.begin(), edges_pool.end());
	edges_pool.erase(unique(edges_pool.begin(), edges_pool.end() ), edges_pool.end());

	for(int j = 0; j < edges_pool.size(); ++j)
	{
		Edge& e = edges_pool[j];
		for (int i = 0; i < 2; i++)
		{
			vertices_pool[e.vertices_id[i]].edges_id.push_back(j);
		}
	}
	
	printf("Loading from %s successfully.\n", f_name);
	printf("Vertex Number = %d\n", vertices_pool.size());
	printf("Triangle Number = %d\n", faces_pool.size());
	printf("Edge Number = %d\n", edges_pool.size());
	fclose(fp);
	return true;
}

void initialize()
{
	// First Update Triangle's matrix
	for(int i = 0; i < faces_pool.size(); ++i)
	{
		Face& f = faces_pool[i];
		f.updateKp();
	}
	// Next Update Point's matrix
	for(int i = 0; i < vertices_pool.size(); ++i)
	{
		Vertex& v = vertices_pool[i];
		v.updateQMatrix();
	}
	// Then Update Edges' err
	for(int i = 0; i < edges_pool.size(); ++i)
	{
		Edge& e = edges_pool[i];
		e.updateNewPoint();
		e.updateErr();
	}
	// Finally Add edges to BST
	for(int i = 0; i < edges_pool.size(); ++i)
	{
		BST.insert(Node(i));
	}
	return;
}

void collapse(int ia, int ib)
{
	Vertex& vb = vertices_pool[ib];
	Vertex& va = vertices_pool[ia];
	vector<int> new_f_index = va.faces_id;
	va.faces_id.clear();
	for(int i = 0; i < vb.faces_id.size(); ++i)
		new_f_index.push_back(vb.faces_id[i]);
	// De-duplicate
	sort(new_f_index.begin(), new_f_index.end());
	new_f_index.erase(unique(new_f_index.begin(), new_f_index.end()), new_f_index.end());
	for(int i = 0; i < new_f_index.size(); ++i)
	{
		Face& f = faces_pool[new_f_index[i]];
		if(!f.enable) continue;
		int j = 0;
		for(; j < 3; ++j)
			if(f.vertices_id[j] == ib)
			{
				f.vertices_id[j] = ia;
				f.confine();
			}		
		if(!f.checkValid())
		{
			f.enable = false;
			continue;
		}
		va.faces_id.push_back(new_f_index[i]);
	}

	vector<int> new_e_index = va.edges_id;
	va.edges_id.clear();
	for(int i = 0; i < vb.edges_id.size(); ++i)
		new_e_index.push_back(vb.edges_id[i]);
	// De-duplicate
	sort(new_e_index.begin(), new_e_index.end());
	new_e_index.erase(unique(new_e_index.begin(), new_e_index.end()), new_e_index.end());
	for(int i = 0; i < new_e_index.size(); ++i)
	{
		Edge& e = edges_pool[new_e_index[i]];
		if(!e.enable) continue;
		int j = 0;
		for(; j < 2; ++j)
			if(e.vertices_id[j] == ib)
			{
				e.vertices_id[j] = ia;
				e.confine();
			}
		
		if(!e.checkValid())
		{
			e.enable = false;
			continue;
		}
		va.edges_id.push_back(new_e_index[i]);
	}

}

void simplify(double ratio)
{
	int number = (1.0 - ratio) / 2.0 * faces_pool.size();
	int turn = 0;
	while (turn++ < number)
	{
		//cout << turn << endl;
		int min_err_edge_id;
		while (1)
		{
			assert(BST.size() > 0);
			min_err_edge_id = BST.begin()->edge_index;
			BST.erase(BST.begin());
			if(edges_pool[min_err_edge_id].enable)
				break;
		}
		Edge& edge = edges_pool[min_err_edge_id];
		int ia = edge.vertices_id[0];
		int ib = edge.vertices_id[1];
		Vector new_pos = edge.v_mid;
		collapse(ia, ib);

		Vertex& va = vertices_pool[ia];
		va.pos = new_pos;
		vertices_pool[ib].enable = false;

		// Now, update Faces of va
		vector<int> update_vertices;
		for(int i = 0; i < va.faces_id.size(); ++i)
		{
			Face& f = faces_pool[va.faces_id[i]];
			assert(f.enable);
			f.updateNormal();
			f.updateKp();
			for(int j = 0; j < 3; ++j)
				update_vertices.push_back(f.vertices_id[j]);
		}
		// De-duplicate
		sort(update_vertices.begin(), update_vertices.end());
		update_vertices.erase(unique(update_vertices.begin(), update_vertices.end()), update_vertices.end());
		// Update Vertices
		vector<int> update_edges;
		for(int i = 0; i < update_vertices.size(); ++i)
		{
			Vertex& v = vertices_pool[update_vertices[i]];
			if(!v.enable)
				continue;
			v.updateQMatrix();
			for(int j = 0; j < v.edges_id.size(); ++j)
				update_edges.push_back(v.edges_id[j]);
		}
		// De-duplicate
		sort(update_edges.begin(), update_edges.end());
		update_edges.erase(unique(update_edges.begin(), update_edges.end()), update_edges.end());
		// Finally, update edges
		for(int i = 0; i < update_edges.size(); ++i)
		{
			Edge& e = edges_pool[update_edges[i]];
			if(!e.enable)
				continue;
			e.updateNewPoint();
			BST.erase(Node(update_edges[i]));
			e.updateErr();
			BST.insert(Node(update_edges[i]));
		}
	}	
}

bool saveOBJFile(const char* f_name)
{
	FILE* fp = fopen(f_name, "w");
	if(fp==NULL)
	{
		printf("Error: Opening %s failed.\n",f_name);
		return false;
	}
	int v_num = 0;
	int f_num = 0;
	map<int, int> dict;
	for(int i = 0; i < vertices_pool.size(); ++i)
	{
		const Vertex& v = vertices_pool[i];
		if(!v.enable)
			continue;
		v_num++;
		dict[i] = v_num;		
		fprintf(fp, "v %lf %lf %lf\n", v.pos[0], v.pos[1], v.pos[2]);
	}
	for(auto f : faces_pool)
	{
		if(!f.enable)
			continue;
		f_num++;
		fprintf(fp, "f %d %d %d\n", dict[f.vertices_id[0]], dict[f.vertices_id[1]], dict[f.vertices_id[2]]);
	}
	fprintf(fp, "# %d %d\n", v_num, f_num);
	return true;
}

int main(int argc, char* argv[])
{
	if (argc != 4)
	{
		cerr << "Input Arguments Invalid.\n";
		cerr << "Usage: MeshSimplify.exe input.obj output.obj 0.3\n";
		return 1;
	}
	string input_file_name(argv[1]);
	string output_file_name(argv[2]);
	double ratio;
	sscanf_s(argv[3], "%lf", &ratio);
	if(!loadOBJFile(input_file_name.c_str()))
	{
		return 1;
	}
	initialize();
	clock_t t1, t2;
	t1 = clock();
	simplify(ratio);
	t2 = clock();
	if(!saveOBJFile(output_file_name.c_str()))
		return 1;
	cout << "Using Time: " << double(t2 - t1) / CLOCKS_PER_SEC << endl;
	return 0;
}




