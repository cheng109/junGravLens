#include <vector> 
#include <climits> 
#include <iostream> 
#include <utility> 
#include <cstdlib> 
#include <ctime> 
using namespace std; 




class Kmeans {
public:
	int k; 
	int n;    // number of data points; 
	vector<int> &group ;  // store the group information; 
	vector<pair<double, double> > centers;  // vector of k centers; 
	vector<double> & xPos; 
	vector<double> & yPos; 

public: 
	Kmeans(int _k, vector<double>& _xPos, vector<double>& _yPos, vector<int>& group); 
	void initiate() ; 
	void updateGroup() ; 
	void updateCenters(); 
	vector<int> countGroup() ; 
}; 









