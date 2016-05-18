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
//#include "gsl/gsl_multimin.h"
//#include <armadillo>
#include <tuple>
#include <map>
//#include <boost/python.hpp>
//#include "gl_crit.h"
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

	string dir 		= argv[1]; 
	string strConf 	= argv[2]; 
	string output 	= argv[3]; 
	/***** prepare*****/
	string confFileName = dir+ strConf; 
	map<string, string> mapConf = parseConfigure(confFileName);


	Image* dataImage = new Image(mapConf["imageFileName"]);
	dataImage->updateFilterImage(mapConf["regionFileName"], stoi(mapConf["usingRegion"]));

	Conf *conf = new Conf(dataImage, mapConf);
	dataImage->updateBackSubtract(conf->back_mean, conf->back_std); 
	dataImage->updateGridPointType();
	dataImage->updateVarList(1, 0.1); // (threshold, var);
	dataImage->invC = dataImage->getVarMatrix();
	
 
	vec d =dataImage->getMatrixD();
	conf->printConfList();

	//cout << "start " << endl; 
	MultModelParam param = MultModelParam(mapConf);
	
	param.printModels();
	param.mix(); 	 // update 'AllMixModels'; 


/*

	for(int i=0; i<6; ++i) {
		string smallMassRegion = "horseshoe_test/reg_" + to_string(i) + ".reg" ; 
		Image* lensImage1 = new Image("horseshoe_test/img_lens.fits"); 
		string imgName = "color_test/f814_cut.fits"; 
		

		Image* dataImage1 = new Image(imgName);
		double 	background = 0.0 ; 
		background = 0.020 ; // f475 
		background = 0.044 ; // f606
		background = 0.0248 ; // f814

		double luminosity =  getMassLuminosity(lensImage1, dataImage1, smallMassRegion,  background) ; 
		cout << "Luminosity for " << imgName << "\t Region: "  << i << "\t "<< luminosity << endl;  
		//delete lensImage1, dataImage1 ; 
	
	} 

*/

	gridSearchVegetti(conf, param,  dataImage, d, dir, output);	

	delete conf; 
	delete dataImage; 
	return 0;
}





