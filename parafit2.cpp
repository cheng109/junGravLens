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
#include <omp.h>

using namespace std;



void mcFit(Conf* conf, MultModelParam param_old, Image* dataImage, string dir, string outputFileName) {
    size_t nLoops(conf->nLoops), ns(dataImage->dataList.size()), nAccept(0), lag(10), writeChkpt(30);
    double weight(0.5), L0;
    //double sstep(0.01*dataImage->d.maxCoeff()), smax(2.*dataImage->d.maxCoeff()), smin(1.*dataImage->d.minCoeff());
    double L(-std::numeric_limits<double>::max());
    double LMax(L);
    double lambdaS = conf->srcRegLevel;
    int jj(0), kk(0);
    bool moveAll(true);

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
    MC mc(model->param, conf->seed);

	ofstream output;
    string out = "mc_chkpt_"+to_string(conf->seed)+".txt";
    if (conf->resume) {
        mc.load(out, model->param, moveAll, L, LMax);
	    output.open(outputFileName, std::ofstream::out | std::ofstream::app);
    } else {
	    output.open(outputFileName);
    }

    for (size_t loop=0; loop<nLoops; ++loop) {
        double step(0.);
        L0 = L;
        if (moveAll) mc.stepPar(model->param);
        else step = mc.stepPar(model->param, jj, kk);
        //cfac2 = mc.stepPar(src, cfac2, iter2);
        model->copyParam(conf, 3);
        //vector<double> penalty = getPenalty2(model, src[3], dataImage, conf, conf->srcRegType);
        vector<double> penalty = getPenalty(model, dataImage, conf, conf->srcRegType);
        L = -penalty[2];

        if (conf->verbose) {
            cout<< loop << ": " << L << " " << L0 << " " << LMax << "  "
                << penalty[0] << " " << penalty[1] << " "
                << model->param.mixAllModels[3][jj].paraList[kk] << " "
                << model->param.mixAllModels[4][jj].paraList[kk] << " "
                << step << endl;
        }
        if (std::isnan(L) || (L<=L0 && mc.random() > exp((L-L0)*weight))) {
            L = L0;
        } else {
            nAccept++;
            if (moveAll) model->copyParam(3,4);
            else model->param.mixAllModels[4][jj].paraList[kk] = model->param.mixAllModels[3][jj].paraList[kk];
            if (L> LMax) {
                if (moveAll) model->copyParam(3,5);
                else model->param.mixAllModels[5][jj].paraList[kk] = model->param.mixAllModels[3][jj].paraList[kk];
                LMax = L;
            }
            if (nAccept % lag == 1) {
                output << model->param.printCurrentModels(4).at(0)
                       << std::scientific << std::setprecision(3) << L
                       << std::setw(10) << std::fixed << loop << endl;
            }
            if (nAccept % writeChkpt == 1) {
                mc.checkPoint(out, model->param, moveAll, L, LMax);
            }
            //for (size_t i=0; i<ns; ++i) src[4](i) = src[3](i);
            //if (L> LMax) {
            //    for (size_t i=0; i<ns; ++i) src[5](i) = src[3](i);
            //    LMax = L;
            //}
        }
    }
    mc.checkPoint(out, model->param, moveAll, L, LMax);
    mc.printIterNum(moveAll);

    output << model->param.printCurrentModels(5).at(0)
           << std::scientific << std::setprecision(3) << LMax
           << std::setw(10) << std::fixed << nLoops << endl;
	output.close();

    model->copyParam(conf, 5);
    vector<double> bestChi = getPenalty(model, dataImage, conf, conf->srcRegType);
    //vector<double> bestChi = getPenalty2(model, src[5], dataImage, conf, conf->srcRegType);
    //for (size_t i=0; i<ns; ++i) model->s[i] = src[5](i);
    cout << "best chi: " << bestChi[0] << " dof: " << ns << " reg: "<< bestChi[1] << " total:" << bestChi[2] << endl;

    writeSrcModResImage(model, dataImage, conf, "mc", dir) ;

    //for (size_t i=0; i<ns; ++i) model->s[i] = sqrt(src[0](i)/src[1](i)-pow(src[2](i)/src[1](i),2));
    //model->writeSrcImage(dir + "img_src_mc_std.fits", conf);


	// Print out the best model :
	cout << "************************\nThe best models : "<< endl;
	cout << model->param.printCurrentModels(5).at(1);
	cout << "************************\n" << endl;
	cout << model->param.printCurrentModels(6).at(1);
	cout << model->param.printCurrentModels(7).at(1);
	cout << model->param.printCurrentModels(8).at(1);
	delete model;

}

Model *model;
#pragma omp threadprivate(model)
void mcFitGW(Conf* conf, MultModelParam param_old, Image* dataImage, string dir, string outputFileName) {
    size_t nLoops(conf->nLoops), lag(5), thin(15), writeChkpt(30), iter(0);
    double lambdaS = conf->srcRegLevel;
    size_t nWalkers(conf->nWalkers);
    auto objective = std::bind(getLogProb, std::placeholders::_1, dataImage, conf);

    #pragma omp parallel
    model = new Model(conf, param_old, lambdaS);

    MC mc(model, conf, objective, nWalkers, outputFileName, iter);

    nLoops += iter;
    for (size_t loop=iter; loop<nLoops; ++loop) {
        #pragma omp parallel for
        for (size_t m=0; m<nWalkers; ++m) {
            mc.stretchMove(model,m);
        }
        if (loop % lag == 0) {
            mc.writeOutput(loop, thin);
            if (loop % writeChkpt == 0) mc.checkPoint(loop);
        }
    }
    mc.writeOutput(nLoops);
    mc.checkPoint(nLoops);
    mc.copyParam(model->param);
    model->copyParam(conf, 5);
    vector<double> bestChi = getPenalty(model, dataImage, conf, conf->srcRegType);
    cout << "best chi: " << bestChi[0] << " dof: " << dataImage->dataList.size() << " reg: "<< bestChi[1] << " total:" << bestChi[2] << endl;

    writeSrcModResImage(model, dataImage, conf, "mcgw_" + to_string(conf->seed), dir) ;

    // Print out the best model :
    cout << "************************\nThe best models : "<< endl;
    cout << model->param.printCurrentModels(5).at(1);
    cout << "************************\n" << endl;
    cout << model->param.printCurrentModels(6).at(1);
    cout << model->param.printCurrentModels(7).at(1);
    cout << model->param.printCurrentModels(8).at(1);
    #pragma omp parallel
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

double getLogProb(Model* model, Image* dataImage, Conf* conf) {

    model->copyParam(conf, 3);
    model->updatePosMapping(dataImage, conf);  // time used: 0.03s;
    model->update_H_zero(conf);
    model->updateLensAndRegularMatrix(dataImage, conf, conf->srcRegType);  // get matrix 'L' and 'RTR'; most time consuming part;
    model->solveSource(&dataImage->invC, &dataImage->d);

    vec &s = model->s;
    vec res = ( model->L * s - dataImage->d);
    vec srcR = s  .transpose() *  model->REG      * s   * model->lambdaS* model->lambdaS  ;

    double lpChi = (res.cwiseProduct(dataImage->invSigma)).squaredNorm()*model->lambdaC* model->lambdaC;
    double lp = srcR[0] + lpChi;
    if (std::isnan(lp)) {
        cout << "Penalty NaN " << lpChi << " " << srcR[0] << endl;
        return -1.0;
    }
    return lp;

}

