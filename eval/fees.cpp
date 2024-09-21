#include "withGraphs.h"
#include <ctime>

#define MAXV 1000000
#define MAXE 20000000
int P;
int Q;

const double feeCustomized = 1.0;
const double feePrescribed = 1.0;

struct PCN {
	int num;
	int edge_count;
	vector<int> adj[MAXV]; 	// this is the vector of edgeIndex, not vector of vertices

	// When representated in theoretical graph objects,
	// each bilateral channel {x,y} is detached into (x,y) and (y,x) for simplicity.
	int sender[MAXE];	
	int receiver[MAXE];	
	int degree[MAXV];
	int countPrescribed[MAXE];
	int countCustomized[MAXE];
	double revenue[MAXV];
	bool isCustomized[MAXE];
	double globalCost;
	
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

void connect2(PCN &pcn, int a, int b, int w, bool isCustomizedChannel) {
 	int ct = ++ pcn.edge_count;
 	pcn.sender[ct] = a;
 	pcn.receiver[ct] = b;
 	pcn.cap[ct] = w;
	pcn.flow[ct] = 0;
	pcn.isCustomized[ct] = isCustomizedChannel;
	pcn.adj[a].push_back(ct);
	pcn.degree[a] ++;
	pcn.countCustomized[ct] = pcn.countPrescribed[ct] = 0;
}

double rand01() {
	return rand() / (double)RAND_MAX;
}

void build_PCN(PCN &pcn, int P, int Q, int weight, double gamma) {
	pcn.num = Q * Q * Q - Q;

	pcn.edge_count = 0;

	pcn.globalCost = 0.0;

	int deglimit = (P + 1) * gamma;
	int degremain = (P + 1)  -  deglimit;

	LPS lps(P, Q);

	for(int i = 0; i < pcn.num; i++) {
		pcn.degree[i] = 0;
		pcn.revenue[i] = 0.0;
	}

	for(int i = 0; i < pcn.num; i++) {
		for(int j = 0; j < deglimit; j++)
			connect2(pcn, i, lps.indexToAdjacency[i][j], weight, false);
	}

	AlbertBarabasi wGraph(pcn.num, degremain - degremain / 2);

	vector< pair<int,int> >::iterator it = wGraph.edges.begin();
	vector< pair<int,int> >::iterator itEd = wGraph.edges.end();
	int bound = degremain * pcn.num / 2;
	int ct = 0;
	for( ; it != itEd; it ++) {
		++ ct;
		if(ct > bound) break;
		connect2(pcn, it->first, it->second, weight, true);
		connect2(pcn, it->second, it->first, weight, true);
	}
}

int qu[MAXE*2];
int d[MAXV];
int path[MAXV]; // array of edgeIndex
int pre[MAXV];


bool find_path(PCN &pcn, int source, int sink) {	//Dinic
	for(int i = 0; i < pcn.num; i++) {
		d[i] = -1;
	}
	qu[0] = source;
	pre[source] = -1;
	d[source] = 0;

	int *l = qu, *r = qu;
	for(; l <= r; l++) {
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
	// when edge weights are equal, BFS already finds the shortest path(s)
	int i = 0, j = sink;
	for(; pre[j] != -1; j = pcn.sender[pre[j]]) {
		path[i++] = pre[j];
		if(pcn.isCustomized[pre[j]]) {
			pcn.countCustomized[pre[j]] ++;
			pcn.revenue[pcn.sender[pre[j]]] += feeCustomized;
			pcn.globalCost += feeCustomized;
		}
		else {
			pcn.countPrescribed[pre[j]] ++;
			pcn.revenue[pcn.sender[pre[j]]] += feePrescribed;
			pcn.globalCost += feePrescribed;
		}
	}
	assert(j == source);
	path[i] = -1;
	return true;
}

PCN pcn_ours;
PCN pcn_typical;


int main() {
	srand(time(NULL));

	double gamma;
	cin >>P >>Q >>gamma;

	if(isQuadraticResidue(P, Q)) {
		cerr <<"(P|Q) should be -1, while P=" <<P <<" and Q=" <<Q <<endl;
		return 0;
	}

	build_PCN(pcn_ours, P, Q, 1, gamma);
	build_PCN(pcn_typical, P, Q, 1, 0.0);

	for(int k = 0; k < 100000; k++) {
		int source = randN() % pcn_ours.num;  // pcn_typical is identical in size
		int sink = randN() % pcn_ours.num;
		if(source == sink) continue;

		find_path(pcn_ours, source, sink);
		find_path(pcn_typical, source, sink);
	}

	vector< pair<int, int> > slice;
	for(int i = 0; i < pcn_typical.num; i++)
		slice.push_back(make_pair(pcn_typical.degree[i], i));
	
	sort(slice.begin(), slice.end(), greater< pair<int,int> >()); // from large to small

	FILE *fout = fopen("feeDistri.csv", "w");

	int ct = 0;

	double A[100], B[100], C[100];

	for(double ratio = .0; ratio < 0.50001; ratio += .02) {
		int index = slice[int(ratio * pcn_typical.num)].second;
		A[ct] = ratio;
		B[ct] = pcn_typical.revenue[index];
		double totalRev = .0;
		for(int sum = 0, i = int(ratio * pcn_typical.num); i < pcn_ours.num; i++) {
			sum += pcn_ours.degree[i];
			if(sum > pcn_typical.degree[index]) {
				sum -= pcn_ours.degree[i];	
				continue;
			}
			totalRev += pcn_ours.revenue[i];
		}
		C[ct] = totalRev;
		ct ++;
	}

	fprintf(fout, "%.2f", A[0]);
	for(int i = 1; i < ct; i++)
		fprintf(fout, ",%.2f", A[i]);
	fprintf(fout, "\n");

	fprintf(fout, "%.0f", B[0]);
	for(int i = 1; i < ct; i++)
		fprintf(fout, ",%.0f", B[i]);
	fprintf(fout, "\n");

	fprintf(fout, "%.0f", C[0]);
	for(int i = 1; i < ct; i++)
		fprintf(fout, ",%.0f", C[i] * pcn_typical.globalCost / pcn_ours.globalCost);
	fprintf(fout, "\n");

	int degree0 = pcn_typical.degree[0];

	cout <<"The global cost of the typical network is: " <<pcn_typical.globalCost <<endl;

	cout <<"The global cost of ours is: " <<pcn_ours.globalCost <<endl;

	fclose(fout);

	return 0;
}
