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
	double lambdaS = 1.0e-5;  

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
		vector<double> penalty = getPenalty(model,  dataImage, conf, "vege") ; 
		
		if(minPenalty > penalty[2]) {
			minPenalty = penalty[2]; 
			minIndex = i; 
		}

		cout << "[" + to_string(i+1) + "/" + to_string(model->param.nComb) + "]\t" << endl; 

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


vector<double> getPenalty(Model* model, Image* dataImage, Conf* conf, string R_type) {
	


	vec d = cV_to_eigenV (&dataImage->dataList); 

	model->updatePosMapping(dataImage, conf);  // time used: 0.03s; 
	model->update_H_zero(conf); 
	model->updateLensAndRegularMatrix(dataImage, conf);  // get matrix 'L' and 'RTR'; most time consuming part; 

	model->solveSource(&dataImage->invC, &dataImage->d, R_type); 
	vec &s = model->s; 
	
	vec res = ( model->L * s - dataImage->d) ; 
	vec chi2 = res.transpose() *  dataImage->invC * res * model->lambdaC* model->lambdaC  ; 
	vec srcR = s  .transpose() *  model->REG      * s   * model->lambdaS* model->lambdaS  ; 



	vector<double> penalty(3); 
	penalty[0] = chi2[0]; 
	penalty[1] = srcR[0]; 
	penalty[2] = chi2[0] + srcR[0]; 

	cout << model->param.parameter[0].critRad << "\t" <<chi2[0] << "\t" << srcR[0] << "\t" << penalty[2] << "\t" << model->H_zero.nonZeros()<< endl; 
	return penalty; 


}


// void gridSearch(Conf* conf, MultModelParam param_old, Image* dataImage, vec d, string dir, string outputFileName) {
// 	Model *model = new Model(conf, param_old, 0.1);
	
// 	//vector<vector<double> > critical;  

// 	vector<int> maxIndex(3,0); 
// 	vector<double> maxObjFunc(3, -1.0); 

// 	int minIndexScatter = 0 ; 
// 	double minScatter = std::numeric_limits<double>::max(); 

// 	ofstream output; 
// 	output.open(outputFileName); 
// 	//MultModelParam newParam (param);  
// 	clock_t	begin = clock(); 

// 	cout << model->param.nComb << endl; 
// 	for(int i=0 ; i< model->param.nComb ; ++i) {
		
// 		model->resetVectors(conf); 
// 		for(int j=0; j<model->param.nLens; ++j) {   // max of j is 3; 
// 			SingleModelParam s; 
// 			s.name = model->param.mixAllModels[i][j].name; 
// 			if(s.name=="PTMASS") {			
// 				s.critRad = model->param.mixAllModels[i][j].paraList[0]; 
// 				s.centerX = model->param.mixAllModels[i][j].paraList[1]; 
// 				s.centerY = model->param.mixAllModels[i][j].paraList[2]; 
// 				model->param.parameter.push_back(s); 
// 			}


// 			if(s.name=="SIE") {
// 				s.critRad = model->param.mixAllModels[i][j].paraList[0]; 
// 				s.centerX = model->param.mixAllModels[i][j].paraList[1]; 
// 				s.centerY = model->param.mixAllModels[i][j].paraList[2]; 
// 				s.e       = model->param.mixAllModels[i][j].paraList[3]; 
// 				s.PA 	  = model->param.mixAllModels[i][j].paraList[4];
// 				s.core 	  = model->param.mixAllModels[i][j].paraList[5];  
// 				model->param.parameter.push_back(s); 
// 			}
// 		}	

// 		vector<double> sBright = dataImage->dataList; 
// 		model->updatePosMapping(dataImage, conf);

// 		double scatter 			= model->getScatterReg(); 
// 		if(scatter < minScatter)  {
// 			minScatter = scatter; 
// 			minIndexScatter = i; 
// 		}

// 		double zerothOrder, gradientOrder, curvatureOrder; 
// 		if (1) {
// 			zerothOrder 		= model->getZerothOrderReg   (conf, sBright);
// 			if(zerothOrder > maxObjFunc[0])  {
// 				maxObjFunc[0] = zerothOrder; 
// 				maxIndex[0] = i; 
// 			}
// 			gradientOrder 	= model->getGradientOrderReg (conf, sBright); 
// 			if(gradientOrder > maxObjFunc[1])  {
// 				maxObjFunc[1] = gradientOrder; 
// 				maxIndex[1] = i; 
// 			}
// 			curvatureOrder 	= model->getCurvatureOrderReg(conf, sBright);
// 			if(curvatureOrder > maxObjFunc[2])  {
// 				maxObjFunc[2] = curvatureOrder; 
// 				maxIndex[2] = i; 
// 			}
// 		}	


// 		if(conf->verbose) {
// 			string pStatus = "[" + to_string(i+1) + "/" + to_string(model->param.nComb) + "]\t" ; 
// 			string resultStatus =  to_string(scatter) + "\t"
// 				+ to_string(zerothOrder) + "\t" 
// 				+ to_string(gradientOrder) + "\t" 
// 				+ to_string(curvatureOrder)  
// 				+ "\n"; 
// 			cout 	<< pStatus << resultStatus ; 
// 			output  << pStatus << model->param.printCurrentModels(i).at(0) << resultStatus ; 		
// 		}

// 		//writeSrcModResImage(model,dataImage,conf, to_string(i), dir) ; 

	
		
// 	}
// 	output.close(); 
	

// 	clock_t end = clock(); 
// 	double elapsed_secs = double(end-begin)/CLOCKS_PER_SEC; 
// 	cout << "Time used: " << elapsed_secs << " seconds. "<< endl; 

// 	// Print out the best model : 
// 	cout << "************************\nThe best models : " << minScatter << endl;
// 	cout << model->param.printCurrentModels(minIndexScatter).at(1);
// 	if(conf->verbose) {
// 		cout << model->param.printCurrentModels(maxIndex[0]).at(1); 
// 		cout << model->param.printCurrentModels(maxIndex[1]).at(1);
// 		cout << model->param.printCurrentModels(maxIndex[2]).at(1); 
// 	}
// 	cout << "************************\n" << endl;
 
// 	delete model;


	
// }



