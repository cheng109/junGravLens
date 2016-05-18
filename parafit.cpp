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

void gridSearchVegetti(Conf* conf, MultModelParam param_old, Image* dataImage, vec d, string dir, string outputFileName) {
	double lambdaS = 10.0; 

	Model *model = new Model(conf, param_old, lambdaS);
	//vector<vector<double> > critical;  

	vector<int> maxIndex(3,0); 
	vector<double> maxObjFunc(3, -1.0); 

	int minIndexScatter = 0 ; 
	double minScatter = std::numeric_limits<double>::max(); 

	ofstream output; 
	output.open(outputFileName); 
	//MultModelParam newParam (param);  
	

	cout << model->param.nComb << endl; 
	for(int i=0 ; i< model->param.nComb ; ++i) {
		
		model->clearVectors (); 
		model->updateReserve(); 
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
		vector<double> sBright = dataImage->dataList; 
		double penalty = getPenalty(&sBright,  model,  dataImage, conf) ; 

		output << model->param.printCurrentModels(i).at(0) << "\t" << penalty << endl; 

	}
	output.close(); 

	// clock_t end = clock(); 
	// double elapsed_secs = double(end-begin)/CLOCKS_PER_SEC; 
	// cout << "Time used: " << elapsed_secs << " seconds. "<< endl; 

	delete model;

}


double getPenalty(vector<double>* sBright, Model* model, Image* dataImage, Conf* conf) {
	


	vec d = cV_to_eigenV (&dataImage->dataList); 

	model->updatePosMapping(dataImage, conf);  // time used: 0.03s; 
	model->updateLensAndRegularMatrix(dataImage, conf);  // get matrix 'L' and 'RTR'; most time consuming part; 
	model->solveSource(&dataImage->invC, d); 
	vec &s = model->s; 
	
	vec res = ( model->L * s - d) ; 
	vec chi2 = res.transpose() *  dataImage->invC * res ; 
	vec srcR = s  .transpose() *  model->HtH      * s   ; 

	double penalty = chi2[0] + model->lambdaS* model->lambdaS * srcR[0]; 

	return penalty; 


}


void gridSearch(Conf* conf, MultModelParam param_old, Image* dataImage, vec d, string dir, string outputFileName) {
	Model *model = new Model(conf, param_old, 0.1);
	
	//vector<vector<double> > critical;  

	vector<int> maxIndex(3,0); 
	vector<double> maxObjFunc(3, -1.0); 

	int minIndexScatter = 0 ; 
	double minScatter = std::numeric_limits<double>::max(); 

	ofstream output; 
	output.open(outputFileName); 
	//MultModelParam newParam (param);  
	clock_t	begin = clock(); 

	cout << model->param.nComb << endl; 
	for(int i=0 ; i< model->param.nComb ; ++i) {
		
		model->clearVectors(); 
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

		vector<double> sBright = dataImage->dataList; 
		model->updatePosMapping(dataImage, conf);





		double scatter 			= model->getScatterReg(); 
		if(scatter < minScatter)  {
			minScatter = scatter; 
			minIndexScatter = i; 
		}

		double zerothOrder, gradientOrder, curvatureOrder; 
		if (1) {
			zerothOrder 		= model->getZerothOrderReg   (conf, sBright);
			if(zerothOrder > maxObjFunc[0])  {
				maxObjFunc[0] = zerothOrder; 
				maxIndex[0] = i; 
			}
			gradientOrder 	= model->getGradientOrderReg (conf, sBright); 
			if(gradientOrder > maxObjFunc[1])  {
				maxObjFunc[1] = gradientOrder; 
				maxIndex[1] = i; 
			}
			curvatureOrder 	= model->getCurvatureOrderReg(conf, sBright);
			if(curvatureOrder > maxObjFunc[2])  {
				maxObjFunc[2] = curvatureOrder; 
				maxIndex[2] = i; 
			}
		}	


		if(conf->verbose) {
			string pStatus = "[" + to_string(i+1) + "/" + to_string(model->param.nComb) + "]\t" ; 
			string resultStatus =  to_string(scatter) + "\t"
				+ to_string(zerothOrder) + "\t" 
				+ to_string(gradientOrder) + "\t" 
				+ to_string(curvatureOrder)  
				+ "\n"; 
			cout 	<< pStatus << resultStatus ; 
			output  << pStatus << model->param.printCurrentModels(i).at(0) << resultStatus ; 		
		}

		if(conf->outputSrcImg) {
			Image* srcImg 	= new Image(model->srcPosXListPixel, model->srcPosYListPixel, &sBright, conf->srcSize[0], conf->srcSize[1], conf->bitpix);		
			string outputSrcName = dir + "img_src_" + to_string(i) +".fits"; 
			if (conf->srcBackground) {
				srcImg -> writeToFile(outputSrcName, conf->back_mean, conf->back_std ) ;
			}
			else 
				srcImg -> writeToFile(outputSrcName) ; 
			delete srcImg; 
		}
		if(conf->outputModImg) {
			//model->updateLensAndRegularMatrix(dataImage, conf);
			//model->updateGradient(dataImage);
			//model->updatePenalty(&dataImage->invC, d);
			Image fullResidualImage = model->getFullResidual(dataImage);
			Image* modelImg = new Image(dataImage->xList, dataImage->yList, &model->mod_img, conf->imgSize[0], conf->imgSize[1], conf->bitpix);
			fullResidualImage.writeToFile(dir + "img_res_" + to_string(i) + ".fits");
			modelImg->writeToFile	(dir + "img_mod_" + to_string(i) +".fits");
			delete modelImg; 
		}
		if(conf->outputCritImg) {
			vector<Image* > curve =  getCritCaustic(conf, &model->param); 
			Image* critImg = curve[0]; 
			Image* causImg = curve[1]; 
			critImg->writeToFile(dir + "img_crit.fits");
			causImg->writeToFile(dir + "img_caus.fits");
			delete critImg; 
			delete causImg; 
		}	
		if(conf->outputLensImg) {
			Image* lensImg = createLensImage(conf, &model->param); 
			lensImg->writeToFile(dir + "img_lens.fits");
			delete lensImg; 
		}	
	}
	output.close(); 
	

	clock_t end = clock(); 
	double elapsed_secs = double(end-begin)/CLOCKS_PER_SEC; 
	cout << "Time used: " << elapsed_secs << " seconds. "<< endl; 

	// Print out the best model : 
	cout << "************************\nThe best models : " << minScatter << endl;
	cout << model->param.printCurrentModels(minIndexScatter).at(1);
	if(conf->verbose) {
		cout << model->param.printCurrentModels(maxIndex[0]).at(1); 
		cout << model->param.printCurrentModels(maxIndex[1]).at(1);
		cout << model->param.printCurrentModels(maxIndex[2]).at(1); 
	}
	cout << "************************\n" << endl;
 
	delete model;


	
}



