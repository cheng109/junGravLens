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
//

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
		}	
		//vector<double> sBright = dataImage->dataList; 
		vector<double> penalty = getPenalty(model,  dataImage, conf) ; 
		
		if(minPenalty > penalty[2]) {
			minPenalty = penalty[2]; 
			minIndex = i; 
		}

		cout << "[" + to_string(i+1) + "/" + to_string(model->param.nComb) + "]\t" ; 
		cout << model->param.parameter[0].critRad << "\t" <<penalty[0] << "\t" << penalty[1] << "\t" << penalty[2] << "\t" <<  endl; 

		writeSrcModResImage(model,dataImage,conf, to_string(i), dir) ; 


		output << model->param.printCurrentModels(i).at(0) << "\t" << penalty[0] <<"\t" <<penalty[1] << "\t" << penalty[2]  << endl; 

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
	
	vec d = cV_to_eigenV (&dataImage->dataList); 

	model->updatePosMapping(dataImage, conf);  // time used: 0.03s; 
	model->update_H_zero(conf); 
	model->updateLensAndRegularMatrix(dataImage, conf);  // get matrix 'L' and 'RTR'; most time consuming part; 
	model->solveSource(&dataImage->invC, &dataImage->d, conf->srcRegType); 
	//model->updateSource(conf) ; // Add noise to the source; 
	vec &s = model->s; 
	vec res = ( model->L * s - dataImage->d) ; 
	vec chi2 = res.transpose() *  dataImage->invC * res * model->lambdaC* model->lambdaC  ; 
	vec srcR = s  .transpose() *  model->REG      * s   * model->lambdaS* model->lambdaS  ; 



	// vector<double> newListX , newListY, srcList; 
	// for(int i=0; i<conf->length; ++i) {
	// 	if (model->srcPosXListPixel[i] <52 and model->srcPosXListPixel[i] > 49
	// 		and model->srcPosYListPixel[i] <60 and model->srcPosYListPixel[i] > 40 ) {
	// 		newListX.push_back(dataImage->xList[i]); 
	// 		newListY.push_back(dataImage->yList[i]); 
	// 		srcList.push_back(1.0); 
	// 		cout << dataImage->xList[i] << "\t" << dataImage->yList[i] << ";\t" ; 

	// 	}

	// }
	//cout << endl; 
	//Image* newImage = new Image(newListX,  newListY, &srcList, conf->imgSize[0], conf->imgSize[1], conf->bitpix); 

	//newImage -> writeToFile ( "check.fits");

	//delete newImage; 

	vector<double> penalty(3); 
	penalty[0] = chi2[0]; 
	penalty[1] = srcR[0]; 
	penalty[2] = chi2[0] + srcR[0]; 

	return penalty; 
}

