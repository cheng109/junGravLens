/*
 * parafit2.cpp
 *
 *  Created on: Apr 22, 2016
 *      Author: En-Hsin Peng
 */


#include "Model.h"
#include "Image.h"
#include <iostream>
#include "commons.h"
#include <fstream>
#include <limits>
#include <iomanip>
#include "parafit.h"
#include "mc.h"

using namespace std;



void mcFit(Conf* conf, MultModelParam param_old, Image* dataImage, string dir, string outputFileName) {
    MC mc(123);
    size_t nLoops(conf->nLoops), ns(dataImage->dataList.size());
    double cfac(1.), weight(0.5), L0;
    //double sstep(0.01*dataImage->d.maxCoeff()), smax(2.*dataImage->d.maxCoeff()), smin(1.*dataImage->d.minCoeff());
    double L(-std::numeric_limits<double>::max());
    double LMax(L);
    double lambdaS = conf->srcRegLevel;
    vector<vector<size_t>> freePar, iter;
    int jj(0), kk(0);

   // vector<vec> src(9, vec(ns));
   // for (size_t i=0; i<9; ++i) src[i].setZero();
   // cout << "max, min: " << smax << "  " << smin << endl;
   // for (size_t i=0; i<ns; ++i) {
   //     double c = (dataImage->dataList[i] > 0) ? dataImage->dataList[i]: 0;
   //     //c = dataImage->dataList[i];
   //     src[3](i) = c;
   //     src[4](i) = c;
   //     src[5](i) = c;
   //     //src[6](i) = smin;
   //     src[7](i) = smax;
   //     src[8](i) = sstep;
   // }
    Model *model = new Model(conf, param_old, lambdaS);
    freePar.resize(model->param.nLens);
    iter.resize(model->param.nLens);
    for(int j=0; j<model->param.nLens; ++j) {
        iter[j].resize(model->param.mixAllModels[0][j].paraList.size(),0);
        for (size_t k=0; k<model->param.mixAllModels[0][j].paraList.size(); ++k) {
            if (model->param.mixAllModels[6][j].paraList[k] < model->param.mixAllModels[7][j].paraList[k]) {
                freePar[j].push_back(k);
            }
        }
    }
	ofstream output;
	output.open(outputFileName);


    for (size_t loop=0; loop<nLoops; ++loop) {
        L0 = L;
        double step = mc.stepPar(model->param, freePar, cfac, iter, jj, kk);
        //cfac2 = mc.stepPar(src, cfac2, iter2);
        model->copyParam(conf, 3);
        //vector<double> penalty = getPenalty2(model, src[3], dataImage, conf, conf->srcRegType);
        vector<double> penalty = getPenalty(model, dataImage, conf, conf->srcRegType);
        L = -penalty[2];

        cout<< loop << ": " << L << " " << L0 << " " << LMax << "  "
            << penalty[0] << " " << penalty[1] << " "
            << model->param.mixAllModels[3][jj].paraList[kk] << " "
            << model->param.mixAllModels[4][jj].paraList[kk] << " "
            << step << endl;
        if (std::isnan(L) || (L<=L0 && mc.random() > exp((L-L0)*weight))) {
            L = L0;
        } else {
            model->param.mixAllModels[4][jj].paraList[kk] = model->param.mixAllModels[3][jj].paraList[kk];
            if (L> LMax) {
                model->param.mixAllModels[5][jj].paraList[kk] = model->param.mixAllModels[3][jj].paraList[kk];
                LMax = L;
            }
            output << model->param.printCurrentModels(4).at(0)
                   << std::scientific << std::setprecision(3) << L
                   << std::setw(10) << std::fixed << loop << endl;
            //for (size_t i=0; i<ns; ++i) src[4](i) = src[3](i);
            //if (L> LMax) {
            //    for (size_t i=0; i<ns; ++i) src[5](i) = src[3](i);
            //    LMax = L;
            //}
        }
    }
    model->copyParam(conf, 5);
    vector<double> bestChi = getPenalty(model, dataImage, conf, conf->srcRegType);
    //vector<double> bestChi = getPenalty2(model, src[5], dataImage, conf, conf->srcRegType);
    //for (size_t i=0; i<ns; ++i) model->s[i] = src[5](i);
    cout << "best chi: " << bestChi[0] << " dof: " << ns << " reg: "<< bestChi[1] << endl;

    cout << "iterations: " << endl;
    for(int j=0; j<model->param.nLens; ++j) {
        cout << j << endl;
        for (auto k: freePar[j]) {
            cout << iter[j][k] << " ";
        }
        cout << endl;
    }

    writeSrcModResImage(model, dataImage, conf, "mc", dir) ;

    //for (size_t i=0; i<ns; ++i) model->s[i] = sqrt(src[0](i)/src[1](i)-pow(src[2](i)/src[1](i),2));
    //model->writeSrcImage(dir + "img_src_mc_std.fits", conf);

	output.close();

	// Print out the best model :
	cout << "************************\nThe best models : "<< endl;
	cout << model->param.printCurrentModels(5).at(1);
	cout << "************************\n" << endl;
	cout << model->param.printCurrentModels(6).at(1);
	cout << model->param.printCurrentModels(7).at(1);
	cout << model->param.printCurrentModels(8).at(1);
	delete model;

}


vector<double> getPenalty2(Model* model, vec &s, Image* dataImage, Conf* conf, string R_type) {

    model->updatePosMapping(dataImage, conf);
    model->update_H_zero(conf);
    model->updateLensAndRegularMatrix(dataImage, conf, R_type);

    vec res = ( model->L * s - dataImage->d) ;
    vec chi2 = res.transpose() *  dataImage->invC * res * model->lambdaC* model->lambdaC;
    vec srcR = s  .transpose() *  model->REG      * s   * model->lambdaS* model->lambdaS;

    vector<double> penalty(3);
    penalty[0] = chi2[0];
    penalty[1] = srcR[0];
    penalty[2] = chi2[0] + srcR[0];

    return penalty;

}

