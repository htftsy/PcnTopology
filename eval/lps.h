#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <map>
#include <iterator>
#include <cassert>
using namespace std;

int signedInverse(int a, int q) {
	// since q is small, it is faster to enumerate directly
	bool opp = false;
	if(a < 0) {
		opp = true;
		a = -a;
	}
	for(int i = 0; i < q; i++)
		if(a * i % q == 1)
			return i * (opp ?-1 :1);
	assert(false);	// should not reach here
}

bool testPrime(int q) {
	assert(q > 1);
	// since q is small, it is faster to enumerate 
	for(int i = 2; i * i <= q; i++)
		if(q % i == 0)
			return false;
	return true;
}

struct matrice {
	int a, b, c, d;
	int q;

	matrice() {}

	matrice(int _a, int _b, int _c, int _d, int _q) {
		a = (_a + 3 * _q * _q) % _q; 
		b = (_b + 3 * _q * _q) % _q; 
		c = (_c + 3 * _q * _q) % _q; 
		d = (_d + 3 * _q * _q) % _q; 
		q = _q;
	}

	void normalize() {
		int mulfac;
		if(a == 0) {
			if(b == 0) {
				if(c == 0) {
					if(d == 0)
						return;
					else
						mulfac = signedInverse(b, q);
				}
				else
					mulfac = signedInverse(c, q);
			}
			else
				mulfac = signedInverse(b, q);
		}
		else
			mulfac = signedInverse(a, q);
		a = a * mulfac % q;
		b = b * mulfac % q;
		c = c * mulfac % q;
		d = d * mulfac % q;
	}

	bool nonZeroDeterminant() {
		return (a * d - b * c) % q != 0;
	}

	void clone(const matrice &t) {
		q = t.q;
		a = (t.a + 3 * q * q) % q;
		b = (t.b + 3 * q * q) % q;
		c = (t.c + 3 * q * q) % q;
		d = (t.d + 3 * q * q) % q;
	}

	void multipleBy(const matrice &t) {
		assert(q == t.q);
		int newa = (a * t.a + b * t.c) % q;
		int newb = (a * t.b + b * t.d) % q;
		int newc = (c * t.a + d * t.c) % q;
		int newd = (c * t.b + d * t.d) % q;
		a = newa; 
		b = newb;
		c = newc;
		d = newd;
		normalize();
	}

	void print() {
		cout <<a <<' ' <<b <<' ' <<c <<' ' <<d <<endl;
	}
};

bool operator < (const matrice &u, const matrice &v) {
	// defined only for the std mapping, 
	// whose internal RB-tree needs a partial order
	if(u.a == v.a)
		if(u.b == v.b)
			if(u.c == v.c)
				return u.d < v.d;
			else
				return u.c < v.c;
		else
			return u.b < v.b;
	else
		return u.a < v.a;
}

vector<matrice> PGLMembers(int q) {
	// returns all members of PGL(2, F_q)
	vector<matrice> vec;
	vec.clear();
	for(int i = 0; i < q; i++)
		for(int j = 0; j < q; j++)
			for(int k = 0; k < q; k++)
				for(int l = 0; l < q; l++) {
					matrice t(i, j, k, l, q);
					if(! t.nonZeroDeterminant())
						continue;
					t.normalize();
					if(t.a == i && t.b == j && t.c == k && t.d == l)
						vec.push_back(t);
				}
	assert((int)vec.size() == q * q * q - q);
	return vec;
}

vector<matrice> generatingSetMembers(int p, int q) {
	int x, y;
	for(x = 0; x < q; x++)
		for(y = 0; y < q; y++)
			if((x * x + y * y + 1) % q == 0)
				goto lp;
	lp: vector<matrice> vec;
	vec.clear();
	for(int a0 = 0; a0 < p; a0++)
		for(int a1 = 1 - p; a1 < p; a1++)
			for(int a2 = 1 - p; a2 < p; a2++)
				for(int a3 = 1 - p; a3 < p; a3++) {
					if(a0 * a0 + a1 * a1 + a2 * a2 + a3 * a3 != p)
						continue;
					if(p % 4 == 1) {
						if((a0 & 1) == 0 || a0 == 0)
							continue;
					}
					else {
						if(!((a0 & 1) == 0 && a0 > 0) && !(a0 == 0 && a1 > 0))
							continue;
					}
					matrice t(a0 + a1 * x + a3 * y, -a1 * y + a2 + a3 * x,
							 -a1 * y - a2 + a3 * x, a0 - a1 * x - a3 * y, q);
					t.normalize();
					vec.push_back(t);
				}
	assert((int)vec.size() == p + 1);
	return vec;
}

bool isQuadraticResidue(int a, int n) {
	// assume that n > 0 and gcd(a,n) = 1
	// since n is small, it is faster to enumerate directly
	for(int i = 1; i < n; i++)
		if(i * i % n == a % n)
			return true;
	return false;
}

struct LPS {
	int p, q;
	int size;

	vector<matrice> matriceNodes;
	vector<matrice> matriceS;

	map<matrice, int> matriceToIndex;

	map<int, vector<int> > indexToAdjacency;

	LPS(int _p, int _q) {
		p = _p;
		q = _q;
		assert(testPrime(p) && testPrime(q));
		assert(p != q && p > 2 && q > 2);
		// To avoid sudden fluctuations of the charts and help to 
		// get a grasp on the tendencies as the graph size grows, 
		// we only work on PGL rather than PSL.
		assert(! isQuadraticResidue(p, q));
		// even q = 97 suffices to model a graph of size roughly 10^6
		assert(p < 100 && q < 100);

		size = q * q * q - q;

		matriceNodes = PGLMembers(q);
		matriceS = generatingSetMembers(p, q);

		assert((int)matriceNodes.size() == size);
		assert((int)matriceS.size() == p + 1);

		for(int i = 0; i < size; i++)
			matriceToIndex[matriceNodes[i]] = i;

		for(int i = 0; i < size; i++) {
			vector<int> vec;
			vec.clear(); 
			for(int j = 0; j <= p; j++) { 
				matrice t;
				t.clone(matriceNodes[i]);
				t.multipleBy(matriceS[j]);

				assert(matriceToIndex.find(t) != matriceToIndex.end());
				int k = matriceToIndex[t];

				vec.push_back(k);
			}
			indexToAdjacency[i] = vec;
		}
	}
};

