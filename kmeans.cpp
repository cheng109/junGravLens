#include <vector> 
#include <climits> 
using namespace std; 




class Kmeans {
private:
	int k; 
	int muX0; 
	int muY0; 
	int n;    // number of data points; 
	vector<int> group ;  // store the group information; 
	vector<pair<int, int> > centers;  // vector of k centers; 

public: 
	Kmeans(int _k, vector<int>& _xPos, vector<int>& _yPos): k(_k), n(_xPos.size())
							group(n) , centers(_k){
		// initite centers; 

		for(int i=0; i<k ; ++i)
			centers[i] = {_xPos[i], _yPos[i]}; 
		
	}
	void updateGroup() {
		for(int i=0; i<n; ++i) {
			int minDist = INT_MAX; 
			for(int j=0; j<k; ++j) {


			}
		}
	}





}; 




int main() {



	return 0; 
}









