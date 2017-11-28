#include <iostream>
#include <string>
#include <vector>
#include <stack>
using namespace std;

// Auxiliary functions.
void print_matrix(const vector<vector<int> > &M);
void find_possible_secondary_structures();
bool forms_pair(char a, char b);
void create(vector<vector<int> > &M, vector<vector<int> > &tracebackM);
void performTraceback(int i, int j, vector<vector<int> > &tracebackM, std::vector<std::pair<int, int> >* result );
string get_format(const std::vector<std::pair<int, int> >& result, int L);

// GLOBALS 
string str = "GGGGGUAUAGCUCAGGGGUAGAGCAUUUGACUGCAGAUCAAGAGGUCCCUGGUUCAAAUCCAGGUGCCCCCU"; 

/* Implementation of Nussinov Maximal Matching Algorithm. */
int main(int argc, char ** argv){
	find_possible_secondary_structures();
	return 0;
}

void create(vector<vector<int> > &M, vector<vector<int> > &tracebackM){
	M.resize(str.length());
	tracebackM.resize(str.length());
	for (size_t i = 0; i < str.length(); ++i)
	{
		M[i].resize(str.length());
		tracebackM[i].resize(str.length());
	}
}

void find_possible_secondary_structures(){
	size_t n = str.size();
	std::vector<std::pair<int, int> > result;
	// str = "-" + str;  
	
	// create nussinov matrix
	vector<vector<int> > M; 
	vector<vector<int> > tracebackM; 
	create(M, tracebackM);

	// initialize matrix
	for (size_t i = 0; i < n; i++){
		M[i][i] = 0;
		M[i][i-1] = 0;
	}

	int max, temp;
	for( int i = str.length()-1; i >= 0; --i){
		for( int j = 0 ; j < str.length(); ++j){
			max = 0;
			if( i >= j ){
				M[i][j] = 0;
				tracebackM[i][j] = -3;
			}
			else{
				temp = M[i+1][j-1] + ( forms_pair(str[j], str[i]) ? 1:0);
				if (temp > max){
					max = temp;
					tracebackM[i][j] = -1 - (forms_pair(str[j], str[i]) ? 1:0);
				}
				for( int k = i; k < j; ++k){
					temp = M[i][k] + M[k+1][j];
					if (temp > max){
						max = temp;
						tracebackM[i][j] = k;
					}
				}
				M[i][j] = max;
			}
		}
	}

	performTraceback(0, str.length()-1, tracebackM, &result);
 	string a = get_format(result, n);
	cout << "done" << endl;
}

void performTraceback(int i, int j, vector<vector<int> > &tracebackM,std::vector<std::pair<int, int> >* result ){
	if(tracebackM[i][j] == -3){
		return;
	}
	else if( tracebackM[i][j] == -2){
		result -> push_back(std::make_pair(i,j));
		performTraceback(i+1, j-1, tracebackM, result);
	}else if (tracebackM[i][j] == -1){
		performTraceback(i+1, j-1, tracebackM, result);
	}else{
		int k = tracebackM[i][j];
		if((i <=k) && (k+1 <= j)){
			performTraceback(i, k, tracebackM, result);
			performTraceback(k+1, j, tracebackM, result);
		}
	}
}	

string get_format(const std::vector<std::pair<int, int> >& v, int L){
	char res[L];
	cout << "" << L << "" << endl;
	for( int index = 0; index < L; index++){
		res[index] = '.';
	}
	cout << res << endl;

	for(int in = 0; in < v.size(); in++){
		res[v[in].first+1] = '(';
		res[v[in].second+1] = ')';
	}
	cout << res << endl;
	return res;
}

bool forms_pair(char a, char b){
	if( (a == 'G' && b == 'C') || (a == 'C' && b == 'G') 
		|| (a == 'A' && b == 'U') || (a == 'U' && b == 'A')
		|| (a == 'G' && b == 'U') || (a == 'U' && b == 'G')){
		return true;
	}
	return false;
}