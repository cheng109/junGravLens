#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <map>
#include <chrono>
#include <ctime>
#include "Image.h"
#include "commons.h"
#include "fitsio.h"
#include "Model.h"
#include "parafit.h"
#include "gsl/gsl_multimin.h"
//#include <armadillo>
#include <tuple>
#include <map>
//#include <boost/python.hpp>
#include "gl_crit.h"
//#include "gl_crit.h"

#include <Eigen/SparseCholesky>
#include <Eigen/Dense>

using namespace std;
//using namespace arma;


int main(int argc, char* argv[]) {


	
	if (argc<4)  { 
		cout << "***************Usage***************" << endl; 
		cout << "3 arguements are required. For example: " << endl; 
		cout << "./junGL horseshoe_test/ conf.txt output.txt" << endl; 
		return 1; 
	}


	// Assign arguments: 

	string dir = argv[1]; 
	string strConf = argv[2]; 
	string output = argv[3]; 
	/***** prepare*****/
	//string dir("pt_test/");
	//string dir("nfw_test/"); 
	//string dir("sersic_test/"); 
	//string dir("sie_test/");
	//("horseshoe_test/");
	//string dir("spemd_test/"); 
	//string dir("blind_test/");
	string confFileName = dir+ strConf; 
	map<string, string> mapConf = parseConfigure(confFileName);


	Image* dataImage = new Image(mapConf["imageFileName"]);
	dataImage->updateFilterImage(mapConf["regionFileName"], stoi(mapConf["usingRegion"]));
	dataImage->updateBackSubtract(stof(mapConf["back_mean"]), stof(mapConf["back_std"])); 
	dataImage->updateGridPointType();
	dataImage->updateVarList(1, 0.1); // (threshold, var);
	dataImage->invC = dataImage->getVarMatrix();
	Conf *conf = new Conf(dataImage, mapConf);
 
	vec d =dataImage->getMatrixD();
	conf->printConfList();

	//cout << "start " << endl; 
	MultModelParam param = MultModelParam(mapConf);
	
	param.printModels();
	param.mix(); 	 // update 'AllMixModels'; 


	gridSearch(conf, param,  dataImage, d, dir, output);	

	delete conf; 
	delete dataImage; 
	return 0;
}





