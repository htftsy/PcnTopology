#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <cmath>
#include <vector>
#include <map>
#include <utility>
#include <iterator>
using namespace std;

struct paras {
	int P;
	int Q;
	double gamma;
	double gamma_origin;

	paras() {}
	paras(int _P, int _Q, double _gamma, double _gamma_origin) 
		:P(_P), Q(_Q), gamma(_gamma), gamma_origin(_gamma_origin) {}
};

int main() {

	vector< pair<int,int> > experiments_PQ;

	experiments_PQ.push_back(make_pair(5, 23));
	experiments_PQ.push_back(make_pair(7, 23));
	experiments_PQ.push_back(make_pair(11, 23));
	experiments_PQ.push_back(make_pair(17, 23));
	experiments_PQ.push_back(make_pair(19, 23));

	experiments_PQ.push_back(make_pair(5, 47));
	experiments_PQ.push_back(make_pair(11, 47));
	experiments_PQ.push_back(make_pair(13, 47));
	experiments_PQ.push_back(make_pair(19, 47));
	experiments_PQ.push_back(make_pair(23, 47)); 
//	experiments_PQ.push_back(make_pair(29, 47));

	vector<paras> experiments;

	for(vector< pair<int,int> >::iterator i = experiments_PQ.begin(); i != experiments_PQ.end(); i++) {
		int m = i->first + 1;
		double slice = 1 / (double) m;
		double epsilon = slice / 2.0;	// to prevent truncations into a smaller integer
		
		for(int j = 0; j <= m; j += 2) // only consider even numbers of open-choice edges
			experiments.push_back(paras(i->first, i->second, slice * j + epsilon, slice * j));
	}

	vector<paras>::iterator st = experiments.begin();
	vector<paras>::iterator ed = experiments.end();
	vector<paras>::iterator it;

	for(it = st; it != ed; it ++) {
			double ave_sum = 0.0;
			double dia_sum = 0.0;
			double tra_sum = 0.0;
			double clu_sum = 0.0;
			double cut_sum = 0.0;
			for(int iterations = 0; iterations < 6; iterations ++){

				FILE *fin = fopen("in.txt", "w");
				fprintf(fin, "%d %d %.2f", it->P, it->Q, it->gamma);

				fclose(fin);
				system("./simu <in.txt>out.txt");
				FILE *fres = fopen("out.txt", "r");
				double ave, dia, tra, clu, cut;
				fscanf(fres, "%lf %lf %lf %lf %lf", &ave, &dia, &tra, &clu, &cut);
				ave_sum += ave;
				dia_sum += dia;
				tra_sum += tra;
				clu_sum += clu;
				cut_sum += cut;
			}
			ave_sum /= 6.0;
			dia_sum /= 6.0;
			tra_sum /= 6.0;
			clu_sum /= 6.0;
			cut_sum /= 6.0;
			cout <<it->P <<',' <<it->Q <<',' <<it->gamma_origin <<','<<ave_sum <<','<<dia_sum <<',' 
				 <<tra_sum <<',' <<clu_sum <<','<<cut_sum <<endl;
	}

}
