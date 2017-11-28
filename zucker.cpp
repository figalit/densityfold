#include <iostream>
#include <string>
#include <vector>

using namespace std;
#define infinity 9999999
#define INF 100
#define a -3.0

// ----------------------------------------------------
// ----------------------------------------------------
// ----------------------------------------------------
// Method for calculating the energies of Hairpin loop.
float eh(char x, char y)
{
	if ((x == 'G' && y == 'C') || (x == 'C' && y == 'G'))
		return -20.4;
	else if ((x == 'A' && y == 'U') || (x == 'U' && y == 'A'))
		return -20.7;
	else if ((x == 'G' && y == 'U') || (x == 'U' && y == 'G'))
		return -20.5;
	else
		return INF;
}

// Method for calculating the energies of bulge or interior loop.
float ebi(char i, char j, char x, char y)
{
	if ((i == 'A' && j == 'U' && x == 'G' && y == 'C') ||
		(i == 'U' && j == 'A' && x == 'C' && y == 'G') ||
		(i == 'A' && j == 'U' && x == 'C' && y == 'G') ||
		(i == 'U' && j == 'A' && x == 'G' && y == 'C'))
		return -8.1;
	else if ((i == 'A' && j == 'U' && x == 'U' && y == 'G') ||
		(i == 'A' && j == 'U' && x == 'G' && y == 'U') ||
		(i == 'U' && j == 'A' && x == 'U' && y == 'G') ||
		(i == 'U' && j == 'A' && x == 'G' && y == 'U'))
		return -8.2;
	else if (( i == 'A' && j == 'U' && x == 'A' && y == 'U') ||
		(i == 'A' && j == 'U' && x == 'U' && y == 'A') ||
		(i == 'U' && j == 'A' && x == 'A' && y == 'U') ||
		(i == 'U' && j == 'A' && x == 'U' && y == 'A'))
		return -8.3;
	else if ((i == 'C' && j == 'G' && x == 'C' && y == 'G') ||
		(i == 'C' && j == 'G' && x == 'G' && y == 'C') ||
		(i == 'G' && j == 'C' && x == 'C' && y == 'G') ||
		(i == 'G' && j == 'C' && x == 'G' && y == 'C'))
		return -8.15;
	else if ((i == 'C' && j == 'G' && x == 'A' && y == 'U') ||
		(i == 'C' && j == 'G' && x == 'U' && y == 'A') ||
		(i == 'G' && j == 'C' && x == 'A' && y == 'U') ||
		(i == 'G' && j == 'C' && x == 'U' && y == 'A'))
		return -8.25;
	else if ((i == 'C' && j == 'G' && x == 'G' && y == 'U') ||
		(i == 'C' && j == 'G' && x == 'U' && y == 'G') ||
		(i == 'G' && j == 'C' && x == 'U' && y == 'G') ||
		(i == 'G' && j == 'C' && x == 'G' && y == 'U'))
		return -8.35;
	else if ((i == 'G' && j == 'U' && x == 'C' && y == 'G') ||
		(i == 'G' && j == 'U' && x == 'G' && y == 'C') ||
		(i == 'U' && j == 'G' && x == 'C' && y == 'G') ||
		(i == 'U' && j == 'G' && x == 'G' && y == 'C'))
		return -8.05;
	else if ((i == 'G' && j == 'U' && x == 'A' && y == 'U') ||
		(i == 'G' && j == 'U' && x == 'U' && y == 'A') ||
		(i == 'U' && j == 'G' && x == 'A' && y == 'U') ||
		(i == 'U' && j == 'G' && x == 'U' && y == 'A'))
		return -8.06;
	else if ((i == 'G' && j == 'U' && x == 'G' && y == 'U') ||
		(i == 'G' && j == 'U' && x == 'U' && y == 'G') ||
		(i == 'U' && j == 'G' && x == 'U' && y == 'G') ||
		(i == 'U' && j == 'G' && x == 'G' && y == 'U'))
		return -8.07;
	else
		return INF;
}

// Method for calculating the energies of Stacking loop.
float es(char x, char y)
{
	if ((x == 'G' && y == 'C') || (x == 'C' && y == 'G'))
		return -12.0;
	else if ((x == 'A' && y == 'U') || (x == 'U' && y == 'A'))
		return -12.1;
	else if ((x == 'G' && y == 'U') || (x == 'U' && y == 'G'))
		return -11.8;
	else
		return INF;
}

// ----------------------------------------------------
// ----------------------------------------------------
// ----------------------------------------------------


/* Declaration of helper functions. */
int find_possible_secondary_structures(string str);

// additional helper functions
void print_matrix(const vector<vector<float> > &A);
float computeV(int i, int j, string str, vector<vector<float> > &V, 
	vector<vector<float> > &W, vector<vector<float> > &tb_W, vector<vector<float> > &tb_V);
void computeW(vector<vector<float> > &W, int L, string str, vector<vector<float> > &V, 
	vector<vector<float> > &tb_W, vector<vector<float> > &tb_V);
// string get_format(vector<Point> v, int L);
void print_matrix(const vector<vector<float> > &A);
void tracebackV(int i, int j, vector<vector<float> > &tb_W, vector<vector<float> > &tb_V, std::vector<std::pair<int, int> >* result);
void tracebackW(int i, int j, vector<vector<float> > &tb_W, vector<vector<float> > &tb_V, std::vector<std::pair<int, int> >* result);
void printStructure(const std::vector<std::pair<int, int> >& result, int L);

/* Simple implementation of sequence alignment between two words using DP algorithm. */
int main(int argc, char **argv){

	string str = "AGACGACAAGGUUGAAUCGCACCCACAGUCUAUGAGUCGGUGACAACAUUACGAAAGGCUGUAAAAUCAAUUAUUCACCACAGGGGGCCCCCGUGUCUAG"; // "GGGAAAUCC";
	find_possible_secondary_structures(str);
	return 0;
}


int find_possible_secondary_structures(string str){

	std::vector<std::pair<int, int> > result;
	size_t n = str.size();
	// str = "-" + str;  
	
	// initialize
	// Note: could also be done line vector<vector<float> > M(n+1, vector<int> (n+1))
	vector<vector<float> > W;
	vector<vector<float> > V;

	W.resize(str.length());
	V.resize(str.length());
	for (size_t i = 0; i < str.length(); ++i)
	{
		V[i].resize(str.length());
		W[i].resize(str.length());
	}

	int m = 3; //4; // enforcing at least how many positions two pase pairs will be far from each other.

	for( int i = 0; i < n; i++){
		for( int j = i+1; j < n ; j++){
			if (i >= j-m) continue; 
			W[i][j] = 0;
			V[i][j] = infinity;
		}
	}
	
	
	vector<vector<float> > tb_W;
	vector<vector<float> > tb_V;

	tb_W.resize(str.length());
	tb_V.resize(str.length());
	for (size_t i = 0; i < str.length(); ++i)
	{
		tb_V[i].resize(str.length());
		tb_W[i].resize(str.length());
	}

	for( int i = 0; i < n; i++){
		for( int j = i+1; j < n ; j++){
			if (i >= j-m) continue; 
			tb_W[i][j] = infinity;
			tb_V[i][j] = infinity;
		}
	}

	computeW(W, n, str, V, tb_W, tb_V);
	// print_matrix(W);

	tracebackW(0, str.length()-1, tb_W, tb_V, &result);
	// print_matrix(tb_W);

	printStructure(result, n);
	cout << "done" << endl;
}

// also add values in traceback 
void computeW(vector<vector<float> > &W, int L, string str, 
	vector<vector<float> > &V, vector<vector<float> > &tb_W, vector<vector<float> > &tb_V){
	float i, j, min;
	
	for( int d = 1; d < L; d++){
		for( i = 0; i < L; i++){
			j = i+d;
			if( j < L){
				min = infinity; // init
				if( W[i+1][j] < min){
					min = W[i+1][j];
					tb_W[i][j] = -1; // case 1
				}
				if( W[i][j-1] < min){
					min = W[i][j-1];
					tb_W[i][j] = -2; // case 2
				}
				if( computeV(i, j, str, V, W, tb_W, tb_V) < min){
					min = computeV(i,j, str, V, W, tb_W, tb_V);
					tb_W[i][j] = -3; // case 3
				}
				for( int k = i+1; k < j; k++){
					int s = W[i][k] + W[k+1][j];
					if( s < min){
						min = s;
						tb_W[i][j] = k; // case 4
					}
				}
				W[i][j] = min;
			} 
		}
	}
}

float computeV(int i, int j, string str, vector<vector<float> > &V,  vector<vector<float> > &W, 
	vector<vector<float> > &tb_W, vector<vector<float> > &tb_V){
	float min = 0;
	if( eh(str[i], str[j]) < min){
		min = eh(str[i], str[j]);
		tb_V[i][j] = -1; // case 1
	}
	if( (es(str[i], str[j]) + V[i+1][j-1]) < min){
		min = es(str[i], str[j]) + V[i+1][j-1];
		tb_V[i][j] = -2; // case 2
	}
	
	for( int k = i+1; k < j-1; k++){
		int temp = W[i+1][k] + W[k+1][j-1] + a;
		if( temp < min){
			min = temp;
			tb_V[i][j] = k; // case 3
		}
	}
	for( int i_prime = i+1; i_prime <= j-2; i_prime++){
		for( int j_prime = i_prime +1; j_prime <= j-1; j_prime++){
			if(i_prime - i + j - j_prime < 2){
				int temp = ebi(str[i], str[j], str[i_prime], str[j_prime]) + V[i_prime][j_prime];
				if( temp < min){
					min = temp;
					tb_V[i][j] = -3; // case 4
				}
			}
		}
	}
	V[i][j] = min;
	return V[i][j];
}
void tracebackW(int i, int j, vector<vector<float> > &tb_W, vector<vector<float> > &tb_V, std::vector<std::pair<int, int> >* result){
	cout << "tracebackW" << endl;
	if(tb_W[i][j] == 0){
		return;
	}
	else if(tb_W[i][j] == -1){ // case 1
		tracebackW(i+1, j, tb_W, tb_V, result);
	}
	else if( tb_W[i][j] == -2){ // case 2
		tracebackW(i, j-1, tb_W, tb_V, result);
	}
	else if(tb_W[i][j] == -3){ // case 3
		tracebackV(i, j, tb_W,tb_V, result);
	}
	else{
		int k = tb_W[i][j];
		if( k >= i+1 && k <= j-1){
			tracebackW(i, k, tb_W, tb_V, result);
			tracebackW(k+1, j, tb_W, tb_V, result);
		}
	}
}

void tracebackV(int i, int j, vector<vector<float> > &tb_W, vector<vector<float> > &tb_V, std::vector<std::pair<int, int> >* result){
	if( tb_V[i][j] == 0){
		return;
	}
	else if(tb_V[i][j] == -1){		
		result -> push_back(std::make_pair(i, j));
	}
	else if(tb_V[i][j] == -2){
		tracebackV(i+1, j-1, tb_W, tb_V, result);
		result -> push_back(std::make_pair(i, j));
	}
	else if(tb_V[i][j] == -3){
		result -> push_back(std::make_pair(i, j));
		tracebackV(i + 1, j - 1, tb_W, tb_V, result);
	}
	else{
		int k = tb_V[i][j];
		if( k >= i+1 && k < j-1){
			// W[i+1][k] + W[k+1][j-1]
			tracebackW(i+1, k, tb_W, tb_V, result);
			tracebackW(k+1, j-1, tb_W, tb_V, result);
		}
	}
}

void printStructure(const std::vector<std::pair<int, int> >& result, int L)
{	
	for (size_t i = 0; i < result.size(); i++){
		// we actually keep the index values in result. 
		cout << "("<< result[i].first + 1 << ", " << result[i].second + 1 << ") ";
	}
	cout << endl;
	
	char res[L];

	for( int index = 0; index < L; index++){
		res[index] = '.';
	}

	for(size_t in = 0; in < result.size(); in++){
		int b = result[in].first;
		int s = result[in].second;
		res[b] = '(';
		res[s] = ')';
	}
	string r = res;
	cout << r.substr(0, L) << endl;

}

void print_matrix(const vector<vector<float> > &A){
	for (int i = 0; i < A.size(); i++){
		for (int j = 0; j < A[i].size()-1; j++){
			if (i <= j){
				cout << A[i][j] << "\t";    
			}else{
				cout << '-' << "\t";
			}
			
		}
		cout << endl;
	}
}

