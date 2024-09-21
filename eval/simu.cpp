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

int qu[MAXE*2];
int d[MAXV];

double average_distance(PCN &pcn, int VERTICES_SAMPLED_FOR_DISTANCE) {
	double sum = 0;
	int sample_num = 0;

	for(int iter = 1; iter < VERTICES_SAMPLED_FOR_DISTANCE; iter++) {
		int source = randN() % pcn.num;
		// BFS is sufficient for outputting shortest paths
		qu[0] = source;
		for(int i = 0; i < pcn.num; i++) 
			d[i] = (1<<30); //simulates infinity
		d[source] = 0;
		int *l = qu, *r = qu;
		for(;l <= r; l++) {
			int cur = *l;
			for(vector<int>::iterator it = pcn.adj[cur].begin(); it != pcn.adj[cur].end(); it++) 
				if(d[pcn.receiver[*it]] > d[cur] + 1) {
					*++r = pcn.receiver[*it];
					d[*r] = d[cur] + 1;
				}
		}
		assert(r-qu <= pcn.num);
		// BFS over
	
		for(int i = 0; i < pcn.num; i++)
			if(i != source) {
				if(d[i] == (1<<30))
					return -1.0;	//returning -1 means infinity (not connected)
				sum += d[i];
				sample_num ++;
			}
	}
	return sum / (double) sample_num;
}

double esti_diameter(PCN &pcn, int ESTI_NUM) {

	double sum = 0;
	for(int est = 0; est < ESTI_NUM; est++) {
		int source = randN() % pcn.num;
		qu[0] = source;
		for(int i = 0; i < pcn.num; i++) 
			d[i] = (1<<30); //simulates infinity
		d[source] = 0;
		int *l = qu, *r = qu;
		for(;l <= r; l++) {
			int cur = *l;
			for(vector<int>::iterator it = pcn.adj[cur].begin(); it != pcn.adj[cur].end(); it++) 
				if(d[pcn.receiver[*it]] > d[cur] + 1) {
					*++r = pcn.receiver[*it];
					d[*r] = d[cur] + 1;
				}
		}
		assert(r-qu <= pcn.num);
		int maxd = 0;
		int posi = -1;
		for(int i = 0; i < pcn.num; i++)
			if(d[i] > maxd) {
				maxd = d[i];
				posi = i;
			}
		if(posi == -1)
			return -1;

		if(maxd > sum)
			sum = maxd;
	}
	return sum >= (1<<30) ?-1 :sum;
}

void add_reserve_edges(PCN &pcn) {
	int existing = pcn.edge_count;
	for(int i = 0; i < existing; i++) {
		int a = pcn.sender[i];
		int b = pcn.receiver[i];

	 	int ct = ++ pcn.edge_count;
	 	pcn.sender[ct] = b;
	 	pcn.receiver[ct] = a;
	 	pcn.cap[ct] = 0;
		pcn.flow[ct] = 0;
		pcn.adj[b].push_back(ct);

		pcn.rev[i] = ct;
		pcn.rev[ct] = i;

		pcn.is_rev[i] = false;
		pcn.is_rev[ct] = true;
	}
}

void reset_cap(PCN &pcn) {
	for(int i = 0; i < pcn.edge_count; i++) {
		if(pcn.is_rev[i])
			pcn.cap[i] = 0;
		else 
			pcn.cap[i] = 1;
	}
}

int path[MAXV]; // array of edgeIndex
int pre[MAXV];


bool find_path(PCN &pcn, int source, int sink) {	//Dinic
	for(int i = 0; i < pcn.num; i++) {
		d[i] = -1;
	}
	qu[0] = source;
	pre[source] = -1;
	d[source] = 0;
	

	#define cost(e) (pcn.is_rev[e] ?-1 :1)

	int *l = qu, *r = qu;
	for(; l <= r; l++) {
	//	cerr<<*l<<endl;
		for(vector<int>::iterator it = pcn.adj[*l].begin(); it != pcn.adj[*l].end(); it++) 
			if(pcn.cap[*it]) {
				int u = pcn.receiver[*it];
				if(d[u] == -1) {
					d[u] = d[*l] + 1;
					*++r = u;
					pre[u] = *it;
				}
			}
	}

	if(d[sink] == -1)
		return false;

	int i = 0, j = sink;
	for(; pre[j] != -1; j = pcn.sender[pre[j]]) {
		path[i++] = pre[j];
	}
	assert(j == source);
	path[i] = -1;
	return true;
}

int maxflow(PCN &pcn, int source, int sink) {
	int sum = 0;
	bool s;

	while(true) {
		s = find_path(pcn, source, sink);
		if(! s)
			break;
		sum += 1;

		for(int i = 0; path[i] != -1; i++) {
			int j = path[i];
			pcn.cap[j] -= 1;	// s=1 anyway
			pcn.cap[pcn.rev[j]] += 1;
		}
	}

	return sum;
}

double expected_mincut(PCN &pcn, int PAIR_SAMPLED) {
	add_reserve_edges(pcn);
	double sum = 0;
	for(int i = 0; i < PAIR_SAMPLED; i++) {
		int source = randN() % pcn.num;
		int sink = randN() % pcn.num;
		while(sink == source) {	// should be different
			sink = randN() % pcn.num;
		}

		reset_cap(pcn);

		int res = maxflow(pcn, source, sink);
		sum += res;	
	}
	
	return sum / (double) PAIR_SAMPLED;
}

double measure_transitivity(PCN &pcn, const set< pair<int,int> > &edges, int WALK_SAMPLED) {
	int edgeSetSize = (int)edges.size();
	int ct = 0;
	for(int sample_num = 0; sample_num < WALK_SAMPLED; sample_num ++) {
		int index1 = randN() % pcn.edge_count;
		int x = pcn.sender[index1];
		int y = pcn.receiver[index1];
		int z = pcn.receiver[pcn.adj[y][randN() % pcn.adj[y].size()]];
		if(edges.find(make_pair(x,z)) != edges.end())
			++ ct;
	}
	return ct / (double) WALK_SAMPLED;
}

double measure_clustering_coeff(PCN &pcn, const set< pair<int,int> > &edges) {
	int ct = 0;
	int demoni = 0;
	for(int i = 0; i < pcn.num; i++) {
		vector<int> adjs = pcn.adj[i];
		for(vector<int>::iterator i = adjs.begin(); i != adjs.end(); i++)
			for(vector<int>::iterator j = i + 1; j != adjs.end(); j++) {
				if(edges.find(make_pair(pcn.receiver[*i], pcn.receiver[*j])) != edges.end())
					++ ct;
				demoni ++;
			}
	}
	return ct / (double) demoni;
}

int Tarjan_count = 0;
int Tarjan_dfn[MAXV];
int Tarjan_low[MAXV];
int Tarjan_fa[MAXV];

void resetTarjanVars(PCN &pcn) {
	Tarjan_count = 0;
	for(int i = 0; i < pcn.num; i++) {
		Tarjan_dfn[i] = Tarjan_low[i] = -1;
		Tarjan_fa[i] = -1;
	}
}

void TarjanSearch(PCN &pcn, int u, int pre, int dep = 0) {
	if(dep > 10000) cout <<dep <<endl;
	++ Tarjan_count;
	Tarjan_dfn[u] = Tarjan_low[u] = Tarjan_count;
	Tarjan_fa[u] = pre;
	vector<int>::iterator st = pcn.adj[u].begin();
	vector<int>::iterator ed = pcn.adj[u].end();
	vector<int>::iterator i;
	for(i = st; i != ed; i++) {
		int target = pcn.receiver[*i];
		if(Tarjan_dfn[target] == -1) {
			TarjanSearch(pcn, target, u, dep + 1);
			if(Tarjan_low[u] > Tarjan_low[target])
				Tarjan_low[u] = Tarjan_low[target];
		}
		else if(pre != target) {
			if(Tarjan_low[u] > Tarjan_dfn[target])
				Tarjan_low[u] = Tarjan_dfn[target];
		}
	}
}

pair<double, double> measure_cut_nodes_cut_edges(PCN &pcn) {
	int cut_node_num = 0;
	int cut_edge_num = 0;

	resetTarjanVars(pcn);
	TarjanSearch(pcn, 0, -1);

	int covers = 0;
	for(int i = 1; i < pcn.num; i++) { // start from 1 instead of 0 
		int tmp = Tarjan_fa[i];
		if(tmp == 0)
			covers ++;
		else {
			if(Tarjan_low[i] >= Tarjan_dfn[tmp])
				++ cut_node_num;
		}
	}
	if(covers > 1)
		++ cut_node_num;

	for(int i = 0; i < pcn.num; i++) { 
		int tmp = Tarjan_fa[i];
		if(tmp != -1 && Tarjan_low[i] > Tarjan_dfn[tmp])
			++ cut_edge_num;
	}

	return make_pair(cut_node_num / (double)pcn.num, cut_edge_num / ((double)pcn.edge_count / 2.0));
}

PCN pcn;

int main() {
	srand(time(NULL));

	double gamma;
	cin >>P >>Q >>gamma;

	if(isQuadraticResidue(P, Q)) {
		cerr <<"(P|Q) should be -1, while P=" <<P <<" and Q=" <<Q <<endl;
		return 0;
	}

	build_PCN(pcn, P, Q, 1, gamma);

	cout <<average_distance(pcn, 30) <<' ';
	cout <<esti_diameter(pcn,30) <<' ';
	
//	pair<double,double> cutMeasures = measure_cut_nodes_cut_edges(pcn);
//	cout <<cutMeasures.first <<' ' <<cutMeasures.second <<' ';
//  Omitted, since there is always no cut node/edge, and the Tarjan DFS often causes stack overflow.

	set< pair<int,int> > edges = pcn.outputEdgeSet();
	cout <<measure_transitivity(pcn, edges, 100000) <<' ';
	cout <<measure_clustering_coeff(pcn, edges) <<' ';

//  always measure mincuts in the final step, since it adds in reversed edges and alters the graph
	cout <<expected_mincut(pcn, 30) <<endl;

	return 0;
}
