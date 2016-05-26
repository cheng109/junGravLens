#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <map>
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
    int opt(0);
    if (argc > 4) opt = atoi(argv[4]);          //0: grid, 1: mcmc

	/***** prepare*****/
	string confFileName = dir+ strConf; 
	map<string, string> mapConf = parseConfigure(confFileName);


	Image* dataImage = new Image(mapConf["imageFileName"]);
	dataImage->updateFilterImage(mapConf["regionFileName"], stoi(mapConf["usingRegion"]));

	Conf *conf = new Conf(dataImage, mapConf);
	dataImage->updateBackSubtract(conf->back_mean, conf->back_std); 
	dataImage->updateGridPointType();
	dataImage->updateVarList(90, conf->back_mean, conf->back_std); // (threshold, var);
	dataImage->invC = dataImage->getVarMatrix();
	
 
	conf->printConfList();

	//cout << "start " << endl; 
	MultModelParam param = MultModelParam(mapConf);
	
	param.printModels();
	param.mix(opt); 	 // update 'AllMixModels'; 



	// Image* maskImg = new Image(dataImage->xList, dataImage->yList, &dataImage->dataList, conf->imgSize[0], conf->imgSize[1], conf->bitpix);
	// maskImg -> writeToFile ("galfit_work/img_mask.fits");
	// delete maskImg; 

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


    if (opt == 0)
        gridSearchVegetti(conf, param, dataImage, dir, output);
    else if (opt == 1)
        mcFit(conf, param, dataImage, dir, output);


	delete conf;
	delete dataImage;
	return 0;
}





