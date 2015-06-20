#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <set>
#include <list>
#include "Matrix.h"
#include "Vector.h"

using namespace std;

const double threshold = 100;


struct Edge
{
	explicit Edge(int x = -1, int y = -1, double e = -1.0) : x(x), y(y), err(e) {}
	int x, y;
	double err;
	bool operator<(const Edge& r) const
	{
		if(x != r.x)
			return x < r.x;
		return y < r.y;
	}
	bool operator==(const Edge& r) const
	{
		return (x == r.x) && (y == r.y);
	}
	bool errIsNull() {return (fabs(err) < 1E-7);}
};

struct PtrCMP
{
	bool operator()(const Edge*& lhs, const Edge*& rhs) {return (*lhs < *rhs);}
};

set<pair<double, Edge> > edges_set;
vector<Edge*> edges;
vector<Vector> vertices;
vector<bool> deleted;
vector<set<Edge*, PtrCMP> > faces;
int faces_number;

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
					vertices.push_back(v);
				}
				else 
				{
					printf("Error: Wrong Number of Values(Should be 3). at Line %d\n",lineNumber);
					return false;
				}	
				break;
			case 'f':				/* face */
				{
					if (faces.empty())
					{
						// First time, initialize.
						int n = vertices.size();
						faces.resize(n);
						deleted.resize(n, false);
					}
					int a, b, c;
					if( fscanf(fp, "%d", &a) ==1 && 
						fscanf(fp, "%d", &b)  ==1 &&
						fscanf(fp, "%d", &c)  ==1 )
					{
						a--;
						b--;
						c--;
						Edge *eee = new Edge(min(b, c), max(b, c));
						edges.push_back(eee);
						faces[a].insert(eee);
						eee = new Edge(min(a, c), max(a, c));
						edges.push_back(eee);
						faces[b].insert(eee);
						eee = new Edge(min(b, a), max(b, a));
						edges.push_back(eee);
						faces[c].insert(eee);
						faces_number++;
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
	printf("Loading from %s successfully.\n", f_name);
	printf("Vertex Number = %d\n", vertices.size());
	printf("Triangle Number = %d\n", faces_number);
	fclose(fp);
	return true;
}

double computeEdgeError(const Edge& e)
{
	// TODO: change v to best point
	Vector v = 0.5 * (vertices[e.x] + vertices[e.y]);
	v[3] = 1.0;
	Matrix kp1;
	for(Edge* vv : faces[e.x])
	{
		int ia, ib, ic;
		Vector va, vb, vc;
		ia = e.x;
		ib = vv->x;
		ic = vv->y;
		va = vertices[ia];
		vb = vertices[ib];
		vc = vertices[ic];
		Vector norm = Vector::normalize(Vector::cross(va - vb, vb - vc));
		norm[3] = -1.0 * (norm[0] * va[0] + norm[1] * va[1] + norm[2] * va[2]);
		Matrix::mul_fast(norm, norm, kp1);
	}
	Matrix kp2;
	for(Edge* vv : faces[e.y])
	{
		int ia, ib, ic;
		Vector va, vb, vc;
		ia = e.y;
		ib = vv->x;
		ic = vv->y;
		va = vertices[ia];
		vb = vertices[ib];
		vc = vertices[ic];
		Vector norm = Vector::normalize(Vector::cross(va - vb, vb - vc));
		norm[3] = -1.0 * (norm[0] * va[0] + norm[1] * va[1] + norm[2] * va[2]);
		Matrix::mul_fast(norm, norm, kp2);
	}
	Matrix kp = kp1 + kp2;
	double err = Matrix::mul(v, kp) * v;
	return err;
}

void initialize()
{
	sort(edges.begin(), edges.end());
	edges.erase(unique(edges.begin(), edges.end()), edges.end());

	for(Edge* e : edges)
	{
		double err = computeEdgeError(*e);
		e->err = err;
		edges_set.insert(make_pair(err, *e));
	}
}

void simplify(double ratio)
{
	int final_number = faces_number * ratio;
	while (faces_number > final_number)
	{
		Edge edge = edges_set.begin()->second;
		edges_set.erase(edges_set.begin());
		// TODO: change the new point
		Vector new_point = 0.5 * (vertices[edge.x] + vertices[edge.y]);
		new_point[3] = 1.0;
		int id = vertices.size();
		vertices.push_back(new_point);

		set<int> changed_id;	// ids remain to be change
		for(auto it = faces[edge.x].begin(); it != faces[edge.x].end(); )
		{
			Edge* ee = *it;
			if(ee->x == edge.y || ee->y == edge.y)
			{
				// Delete this Triangle
				if(edge.x < ee->x)
					faces_number--;		// Delete only once
				faces[edge.x].erase(it++);
				continue;
			}
			int ia, ib, ic;
			ia = ee->x;
			ib = ee->y;
			ic = edge.x;
			changed_id.insert(ia);
			changed_id.insert(ib);
			// First remove this Triangle in three vertices
			faces[edge.x].erase(it++);
			Edge eee(min(ib, ic), max(ib, ic));
			faces[ia].erase(&eee);
			eee = Edge(min(ia, ic), max(ia, ic));
			faces[ib].erase(&eee);
			// And then add them back
			Edge* tmp = new Edge(min(ib, id), max(ib, id));
			faces[ia].insert(tmp);
			//edges.push_back(tmp);
			tmp = new Edge(min(ia, id), max(ia, id));
			//edges.push_back(tmp);
			faces[ib].insert(tmp);
		}
		for(auto it = faces[edge.y].begin(); it != faces[edge.y].end(); )
		{
			Edge* ee = *it;
			if(ee->x == edge.x)
			{
				// Delete this Triangle
				faces_number--;
				faces[edge.y].erase(it++);
				continue;
			}
			int ia, ib, ic;
			ia = ee->x;
			ib = ee->y;
			ic = edge.y;
			changed_id.insert(ia);
			changed_id.insert(ib);
			// First remove this Triangle in three vertices
			faces[edge.y].erase(it++);

			Edge eee(min(ib, ic), max(ib, ic));
			faces[ia].erase(&eee);
			eee = Edge(min(ia, ic), max(ia, ic));
			faces[ib].erase(&eee);
			// And then add them back
			Edge* tmp = new Edge(min(ib, id), max(ib, id));
			faces[ia].insert(tmp);
			//edges.push_back(tmp);
			tmp = new Edge(min(ia, id), max(ia, id));
			//edges.push_back(tmp);
			faces[ib].insert(tmp);
		}
		// And now, update BBST
		// First get all affected edges
		set<Edge> changed_edges;
		for(auto i : changed_id)
		{
			for (auto ee : faces[i])
			{
				int ia = ee->x;
				int ib = ee->y;
				changed_edges.insert(Edge(min(i, ia), max(i, ia), ee->err));
				changed_edges.insert(Edge(min(i, ib), max(i, ib), ee->err));
			}
		}
		// Then, update BBST
		for(Edge* ee :changed_edges)
		{
			// First remove edges in BBST
			if(ee)
			
		}


		deleted[edge.x] = true;
		deleted[edge.y] = true;		
	}
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
	sscanf_s(argv[3], "%llf", &ratio);
	loadOBJFile(input_file_name.c_str());
	initialize();
	simplify(ratio);

	set<int> s;
	s.insert(2);
	s.insert(3);
	s.insert(1);
	s.insert(10);

	s.erase(10);
	for(auto i : s)
		cout << i << endl;

	return 0;
}

