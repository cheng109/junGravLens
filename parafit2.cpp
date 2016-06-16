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
    size_t nLoops(1000), iter1(0), iter2(0), ns(dataImage->dataList.size());
    double cfac(1.), cfac2(1.), weight(0.5), L0;
    double sstep(0.01*dataImage->d.maxCoeff()), smax(2.*dataImage->d.maxCoeff()), smin(1.*dataImage->d.minCoeff());
    double L(-std::numeric_limits<double>::max());
    double LMax(L);
    double lambdaS = 1e-5;
    vector<vec> src(9, vec(ns));

    for (size_t i=0; i<9; ++i) src[i].setZero();
    cout << "max, min: " << smax << "  " << smin << endl;
    for (size_t i=0; i<ns; ++i) {
        double c = (dataImage->dataList[i] > 0) ? dataImage->dataList[i]: 0;
        //c = dataImage->dataList[i];
        src[3](i) = c;
        src[4](i) = c;
        src[5](i) = c;
        //src[6](i) = smin;
        src[7](i) = smax;
        src[8](i) = sstep;
    }
    Model *model = new Model(conf, param_old, lambdaS);
	ofstream output;
	output.open(outputFileName);

    for (size_t loop=0; loop<nLoops; ++loop) {
        L0 = L;
        if (loop & 1) {
            cfac = mc.stepPar(model->param, cfac, iter1);
            iter1++;
        } else {
            cfac2 = mc.stepPar(src, cfac2, iter2);
            iter2++;
        }
        model->copyParam(conf, 3);
        vector<double> penalty = getPenalty2(model, src[3], dataImage, conf, "vege");
        L = -penalty[2];

        if (loop & 1) cout<< loop << ": " << L << " " << L0 << " " << LMax << "  "
                          << penalty[0] << " " << penalty[1] << " "
                          << model->param.mixAllModels[3][0].paraList[0] << " "
                          << model->param.mixAllModels[4][0].paraList[0] << " "
                          << model->param.mixAllModels[8][0].paraList[0] << endl;
        if (std::isnan(L) || (L<=L0 && mc.random() > exp((L-L0)*weight))) {
            L = L0;
        } else {
            if (loop & 1) {
                model->copyParam(3, 4);
                if (L> LMax) {
                    model->copyParam(3, 5);
                    LMax = L;
                }
                output << model->param.printCurrentModels(4).at(0)
                       << std::scientific << std::setprecision(3) << L << endl;
            } else {
                for (size_t i=0; i<ns; ++i) src[4](i) = src[3](i);
                if (L> LMax) {
                    for (size_t i=0; i<ns; ++i) src[5](i) = src[3](i);
                    LMax = L;
                }
            }
        }
    }
    model->copyParam(conf, 5);
    vector<double> bestChi = getPenalty2(model, src[5], dataImage, conf, "vege");
    for (size_t i=0; i<ns; ++i) model->s[i] = src[5](i);
    cout << "best chi: " << bestChi[0] << " dof: " << src[0].size() << " reg: "<< bestChi[1] << endl;
    writeSrcModResImage(model, dataImage, conf, "mc", dir) ;

    for (size_t i=0; i<ns; ++i) model->s[i] = sqrt(src[0](i)/src[1](i)-pow(src[2](i)/src[1](i),2));
    model->writeSrcImage(dir + "img_src_mc_std.fits", conf);

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

