#include "withGraphs.h"
#include <ctime>

#define MAXV 1000000
#define MAXE 20000000
int P;
int Q;

struct PCN {
	int num;
	int edge_count;
	vector<int> adj[MAXV]; 	// this is the vector of edgeIndex, not vector of vertices

	// When representated in theoretical graph objects,
	// each bilateral channel {x,y} is detached into (x,y) and (y,x) for simplicity.
	int sender[MAXE];	
	int receiver[MAXE];	
	
	// variables for network flows
	int rev[MAXE];  // the reverse edge, for maxflow only
	bool is_rev[MAXE];
	int cap[MAXE];
	int flow[MAXE];

	set< pair<int,int> > outputEdgeSet() {
		set< pair<int,int> > s;
		s.clear();
		for(int i = 0; i < edge_count; i++)
			s.insert(make_pair(sender[i], receiver[i]));
		return s;
	}
};

void connect2(PCN &pcn, int a, int b, int w) {
 	int ct = ++ pcn.edge_count;
 	pcn.sender[ct] = a;
 	pcn.receiver[ct] = b;
 	pcn.cap[ct] = w;
	pcn.flow[ct] = 0;
	pcn.adj[a].push_back(ct);
	// cout <<a <<' '<<b<<endl;
}

double rand01() {
	return rand() / (double)RAND_MAX;
}

void build_PCN(PCN &pcn, int P, int Q, int weight, double gamma) {
	pcn.num = Q * Q * Q - Q;

	pcn.edge_count = 0;

	int deglimit = (P + 1) * gamma;
	int degremain = (P + 1)  -  deglimit;

	LPS lps(P, Q);

	for(int i = 0; i < pcn.num; i++) {
		for(int j = 0; j < deglimit; j++)
			connect2(pcn, i, lps.indexToAdjacency[i][j], weight);
	}

	AlbertBarabasi wGraph(pcn.num, degremain - degremain / 2);

	vector< pair<int,int> >::iterator it = wGraph.edges.begin();
	vector< pair<int,int> >::iterator itEd = wGraph.edges.end();
	int bound = degremain * pcn.num / 2;
	int ct = 0;
	for( ; it != itEd; it ++) {
		++ ct;
		if(ct > bound) break;
		connect2(pcn, it->first, it->second, weight);
		connect2(pcn, it->second, it->first, weight);
	}
}

PCN pcn;

int main() {
	build_PCN(5, 23, );
}

