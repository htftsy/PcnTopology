#include "lps.h"
#include <set>
#include <climits>

const double watts_parameter = 0.75;

inline int randN() {	// returns a random unsigned int of 30 bits
	return ((rand()&0x7fff)<<15) + (rand()&0x7fff);
}

struct Watts {
	int N; 
	int K1; // K = 2 * K1
	vector< pair<int,int> > edges;
	set< pair<int,int> > edgeInSet;

	bool found_collision(int x, int y) {
		if(x == y)
			return true;
		return edgeInSet.find(make_pair(x, y)) != edgeInSet.end() 
			|| edgeInSet.find(make_pair(y, x)) != edgeInSet.end();
	}

	Watts(int _N, int _K1) :N(_N), K1(_K1) {
		edges.clear();
		edgeInSet.clear();

		for(int i = 0; i < N; i++) {
			for(int j = 1; j <= K1; j++)
				edges.push_back(make_pair(i, (i + j) % N));
		}

		for(vector< pair<int,int> >::iterator i = edges.begin(); i != edges.end(); i++)
			edgeInSet.insert(*i);

		for(int i = 0; i < N; i++) {
			for(int j = 0; j < K1; j++) {
				int edgeIndex = i * K1 + j;
				if(rand() < watts_parameter * (double) RAND_MAX) {
					int target = randN() % N;
					while(found_collision(i, target)) {
						target = randN() % N;
					}
					edgeInSet.erase(edges[edgeIndex]);
					edges[edgeIndex].second = target;
					edgeInSet.insert(edges[edgeIndex]);
				}
			}
		}
	}

	void printAllEdgesInCErr() {
		for(vector< pair<int,int> >::iterator i = edges.begin(); i != edges.end(); i++)
			cerr <<"edge "<<i->first <<' ' <<i->second <<endl;
	}
};

struct AlbertBarabasi {
	int N;
	int K1;
	vector< pair<int,int> > edges;
	vector<int> deg;
	int totalDeg;

	void addEdge(int x, int y) {
		edges.push_back(make_pair(x, y));
		deg[x] = deg[x] + 1;
		deg[y] = deg[y] + 1;
		totalDeg += 2;
	}

	int randomTar(int st, int ed, int totalDeg) {
		int r = randN();
		for(int ct = 0, i = st; i < ed; i++) {
			ct += deg[i];
			if(ct / (double)totalDeg > r / double(1<<30))
				return i;
		}
		return ed - 1;
	}

	AlbertBarabasi(int _N, int _K1) :N(_N), K1(_K1) {
		deg.clear();
		deg.push_back(2);  // deg[0]=2, an ideal (not real) self-loop for self-containness
		totalDeg = 2;
		for(int i = 1; i < N; i++)
			deg.push_back(0);
		for(int i = 1; i < N; i++) {
			for(int j = 0; j < K1; j++) {
				int target = randomTar(0, i, totalDeg);
				addEdge(i, target);
			}
		}
	}

	void printAllEdgesInCErr() {
		for(vector< pair<int,int> >::iterator i = edges.begin(); i != edges.end(); i++)
			cerr <<"edge "<<i->first <<' ' <<i->second <<endl;
		for(int i = 0; i < N; i++) 
			cerr <<i <<" has degree of "<<deg[i] <<endl;
	}
};
