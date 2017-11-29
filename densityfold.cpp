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
float eH(char x, char y)
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
float eBI(char i, char j, char x, char y)
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
float eS(char x, char y)
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

float eM();
float eDA();

// ----------------------------------------------------
// ----------------------------------------------------
// ----------------------------------------------------

// DP tables to be computed. 
vector<float> ED;
vector<float> E;
vector<vector<float> > EDs;
vector<vector<float> > Es;
vector<vector<float> > EDbi;
vector<vector<float> > Ebi;
vector<vector<float> > EDm;
vector<vector<float> > Em;

vector<float> tracebackED;
vector<float> tracebackE;
vector<vector<float> > tracebackEDs;
vector<vector<float> > tracebackEs;
vector<vector<float> > tracebackEDbi;
vector<vector<float> > tracebackEbi; // no need, only one option for traceback
vector<vector<float> > tracebackEDm;
vector<vector<float> > tracebackEm;

// Variables
std::vector<std::pair<int, int> > result;
string str = "GCACGACG"; // "GGGAAAUCC";

// Method declarations.
void findSecondaryStructure(string s);
void initializeTables(int n);
void computeED_j(int n);
float computeEbi(int i, int j);
float computeEm(int i, int j);
float computeEDs(int i, int j);
float computeEs(int i, int j);
float computeEDbi(int i, int j);
float computeEDm(int i, int j);
void print_matrix(const vector<vector<float> > &A);
void printStructure(const std::vector<std::pair<int, int> >& result, int L);
/* 	Initial implementation of density fold 
*	algorithm for RNA structure prediction. 
*/
int main(int argc, char **argv){
	findSecondaryStructure(str);
	return 0;
}

void findSecondaryStructure(string str){
	int n = str.length();
	std::vector<std::pair<int, int> > result;

	initializeTables(n);
 	computeED_j(n);
 	print_matrix(EDs);
 	printStructure(result, n);
}

void computeED_j(int n){
	float min;
	for( int j = 0; j < n; j++){
		min = 0;
		// case 1
		if( ED[j-1] < min){
			min = ED[j-1];
			tracebackED[j] = -1; 
		} 
		// case 2	
		for( int i = 1; i <= j-1; i++){
			float EDs_i_j = computeEDs(i, j);
			if( (ED[i-1] + EDs_i_j) < min){
				float case2_min = EDs_i_j + ED[i-1];
				tracebackED[j] = i; 
			}
		}
	
	}
}

float computeEDs(int i, int j){
	float min;
	// case 1
	min = infinity;
	// case 2
	if (eH(str[i], str[j]) < min){
		min = eH(str[i], str[j]);
		tracebackEDs[i][j] = -2;
	}
	// case 3
	float from_Es = computeEs(i+1, j-1);
	float density = (2*((eS(str[i], str[j]) + from_Es)/(j-i+1)) ) + EDs[i+1][j-1];
	if( density < min){
		min = density;
		tracebackEDs[i][j] = -3; 
	}
	// case 4
	if( computeEDbi(i, j) < min){
		min = computeEDbi(i, j);
		tracebackEDs[i][j] = -4;
	}
	// case 5
	if( computeEDm(i, j) < min){
		min = computeEDm(i, j);
		tracebackEDs[i][j] = -5;
	}
	EDs[i][j] = min;
	return min;
}

float computeEDbi(int i, int j){
	float min = 0;
	for( int i_prime = i+1; i_prime <= j-2; i_prime++){
		for( int j_prime = i_prime +1; j_prime <= j-1; j_prime++){
			if(i_prime - i + j - j_prime < 2){
				int temp = eBI(str[i], str[j], str[i_prime], str[j_prime]) + EDs[i_prime][j_prime];
				if( temp < min){
					min = temp;
					tracebackEDs[i][j] = -3; // case 4
				}
			}
		}
	}
	EDbi[i][j] = min;
	return min;
}

// TODO
float computeEDm(int i, int j){
	return -3.0;
}

float computeEs(int i, int j){
	float min;
	// case 1
	min = infinity;
	// case 2
	if(eH(str[i], str[j]) < min){
		min = eH(str[i], str[j]);
		tracebackEs[i][j] = -2;
	}
	// case 3
	if (i < j)
	if( (Es[i+1][j-1] + eS(str[i], str[j])) < min){
		min = Es[i+1][j-1] + eS(str[i], str[j]);
		tracebackEs[i][j] = -3;
	}
	// case 4
	if( computeEbi(i, j) < min){
		min = computeEbi(i, j);
		tracebackEs[i][j] = -4;
	}
	// case 5
	if( computeEm(i, j) < min){
		min = computeEm(i, j);
		tracebackEs[i][j] = -5;
	}
	Es[i][j] = min;
	return min;
}

float computeEbi(int i, int j){
	float min = 0;
	for( int i_prime = i+1; i_prime <= j-2; i_prime++){
		for( int j_prime = i_prime +1; j_prime <= j-1; j_prime++){
			// TODO MAYBE EXTRA CONDITION HERE
			int temp = eBI(str[i], str[j], str[i_prime], str[j_prime]) + Es[i_prime][j_prime];
			if( temp < min){
				min = temp;
			}
		}
	}
	return min;
}

// TODO
float computeEm(int i, int j){
	return -3.0;
}


void initializeTables(int n){
	
	ED.resize(n);
	E.resize(n);
	EDs.resize(n);
	Es.resize(n);
	EDbi.resize(n);
	Ebi.resize(n);
	EDm.resize(n);
	Em.resize(n);

	tracebackED.resize(n);
	tracebackE.resize(n);
	tracebackEDs.resize(n);
	tracebackEs.resize(n);
	tracebackEDbi.resize(n);
	tracebackEbi.resize(n);
	tracebackEDm.resize(n);
	tracebackEm.resize(n);
	
	for (int i = 0; i < n; ++i)
	{
		EDs[i].resize(n);
		Es[i].resize(n);
		EDbi[i].resize(n);
		Ebi[i].resize(n);
		EDm[i].resize(n);
		Em[i].resize(n);

		tracebackEDs[i].resize(n);
		tracebackEs[i].resize(n);
		tracebackEDbi[i].resize(n);
		tracebackEbi[i].resize(n);
		tracebackEDm[i].resize(n);
		tracebackEm[i].resize(n);
	}
}

void printStructure(const std::vector<std::pair<int, int> >& result, int L){	
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