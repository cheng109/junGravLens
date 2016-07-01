/*
 * parafit.cpp
 *
 *  Created on: Dec 24, 2015
 *      Author: cheng109
 */


//#include "gsl/gsl_multimin.h"
#include "Model.h"
#include "Image.h"
#include <iostream>
#include "commons.h"
#include <fstream>
#include <limits>
#include <ctime> 
#include "parafit.h"

using namespace std;

void gridSearchVegetti(Conf* conf, MultModelParam param_old, Image* dataImage, string dir, string outputFileName) {
	double lambdaS = conf->srcRegLevel;  

	Model *model = new Model(conf, param_old, lambdaS);
		
	int minIndex = 0; 
	double minPenalty = std::numeric_limits<double>::max(); 

	ofstream output; 
	output.open(outputFileName); 
	

	cout << model->param.nComb << endl; 
	for(int i=0 ; i< model->param.nComb ; ++i) {
		
		model->resetVectors (conf); 
		//model->updateReserve(conf); 

		for(int j=0; j<model->param.nLens; ++j) {   // max of j is 3; 
			SingleModelParam s; 
			s.name = model->param.mixAllModels[i][j].name; 	
			if(s.name=="PTMASS") {			
				s.critRad = model->param.mixAllModels[i][j].paraList[0]; 
				s.centerX = model->param.mixAllModels[i][j].paraList[1]; 
				s.centerY = model->param.mixAllModels[i][j].paraList[2]; 
				model->param.parameter.push_back(s); 
			}
			if(s.name=="SIE") {
				s.critRad = model->param.mixAllModels[i][j].paraList[0]; 
				s.centerX = model->param.mixAllModels[i][j].paraList[1]; 
				s.centerY = model->param.mixAllModels[i][j].paraList[2]; 
				s.e       = model->param.mixAllModels[i][j].paraList[3]; 
				s.PA 	  = model->param.mixAllModels[i][j].paraList[4];
				s.core 	  = model->param.mixAllModels[i][j].paraList[5];  
				model->param.parameter.push_back(s); 
			}

			if(s.name=="NFW") {
				s.massScale = model->param.mixAllModels[i][j].paraList[0]; 
				s.centerX   = model->param.mixAllModels[i][j].paraList[1]; 
				s.centerY   = model->param.mixAllModels[i][j].paraList[2]; 
				s.e         = model->param.mixAllModels[i][j].paraList[3]; 
				s.PA 	    = model->param.mixAllModels[i][j].paraList[4];
				s.radScale  = model->param.mixAllModels[i][j].paraList[5];  
				model->param.parameter.push_back(s); 
			}

			if(s.name=="SPEMD") {
				s.critRad = model->param.mixAllModels[i][j].paraList[0]; 
				s.centerX = model->param.mixAllModels[i][j].paraList[1]; 
				s.centerY = model->param.mixAllModels[i][j].paraList[2]; 
				s.e       = model->param.mixAllModels[i][j].paraList[3]; 
				s.PA 	  = model->param.mixAllModels[i][j].paraList[4];
				s.core 	  = model->param.mixAllModels[i][j].paraList[5]; 
				s.power	  = model->param.mixAllModels[i][j].paraList[6]; 
				model->param.parameter.push_back(s); 
			}

		}	
		vector<double> penalty = getPenalty(model,  dataImage, conf) ; 
		
		if(minPenalty > penalty[2]) {
			minPenalty = penalty[2]; 
			minIndex = i; 
		}

		cout << "[" + to_string(i+1) + "/" + to_string(model->param.nComb) + "] " ; 
		//cout <<"\t" << model->param.parameter[0].critRad << "\t" <<penalty[0] << "\t" << penalty[1] << "\t" << penalty[2] << "\t"; 
		
		writeSrcModResImage(model,dataImage,conf, to_string(i), dir) ; 
		//output << model->param.printCurrentModels(i).at(0) << "\t" << penalty[0] <<"\t" <<penalty[1] << "\t" << penalty[2]  << endl; 
		cout << endl; 
	}

	cout << "************************\nThe best models : " << minPenalty << endl;
	cout << model->param.printCurrentModels(minIndex).at(1);
	cout << "************************\n" << endl;	

	output.close(); 

	// clock_t end = clock(); 
	// double elapsed_secs = double(end-begin)/CLOCKS_PER_SEC; 
	// cout << "Time used: " << elapsed_secs << " seconds. "<< endl; 

	delete model;

}


vector<double> getPenalty(Model* model, Image* dataImage, Conf* conf) {
	vector<double> penalty(3); 
	vec d = cV_to_eigenV (&dataImage->dataList); 

	model->updatePosMapping(dataImage, conf);  // time used: 0.03s; 
	if(1) {
		model->update_H_zero(conf); 
		model->updateLensAndRegularMatrix(dataImage, conf);  // get matrix 'L' and 'RTR'; most time consuming part; 
		model->solveSource(&dataImage->invC, &dataImage->d, conf->srcRegType); 

		// modify source ; srcPosXListPixel ; 
		double x0 = 54; 
		double y0 = 57; 
		double  a = 1.5; 
		double  b = y0-x0*a; 
		for(int i=0; i<model->srcPosXListPixel.size(); ++i) {
			int x = model->srcPosXListPixel[i]; 
			int y = model->srcPosYListPixel[i]; 
			if( y - (a*x + b) <= 0 ) {
				//model->s[i] = 0; 
			}
		}


		vec &s = model->s; 
		vec res = ( model->L * s - dataImage->d) ; 
		vec chi2 = res.transpose() *  dataImage->invC * res * model->lambdaC* model->lambdaC  ; 
		vec srcR = s  .transpose() *  model->REG      * s   * model->lambdaS* model->lambdaS  ; 



		penalty[0] = chi2[0]; 
		penalty[1] = srcR[0]; 
		penalty[2] = chi2[0] + srcR[0]; 
	}
	//penalty[2] = model->getScatterReg() ; 

	return penalty; 
}

