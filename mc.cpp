/*
 * mc.cpp
 *
 *  Created on: Apr 23, 2016
 *      Author: En-Hsin Peng
 */

#include "mc.h"
#include <iostream>

MC::MC(unsigned seed) {
    rng_engine.seed(seed);
    rng = std::bind(std::uniform_real_distribution<double>(0.,1.), std::ref(rng_engine));
    normal = std::bind(std::normal_distribution<double>(0.,1.), std::ref(rng_engine));
    //makeCgauss();
}

double MC::random() {
    return rng();
}

double MC::gauss(double sigma, double mu, double x) {
    double dx = x- mu;
    return exp(-dx*dx/(2*sigma*sigma))/(sigma*sqrt(2*M_PI));
}

void MC::makeCgauss() {
    cgArr.resize(100,0.);
    cgArr[0] = gauss(1., 0.,  -5.)*0.1;
    for (size_t i=1; i<100; ++i) {
        cgArr[i] = cgArr[i-1] + gauss(1., 0.,  -5. + i*0.1)*0.1;
    }
}

double MC::cgauss() {
    //double xv = 0.;
    //double r = rng();
    //for (size_t i=0; i<100; ++i) {
    //    if (cgArr[i]>r) {
    //        xv = (-5. + i*0.1) - 0.1*rng();
    //        break;
    //    }
    //}
    //return xv;
    return normal();
}


double MC::stepPar(MultModelParam &param, vector<vector<size_t>> &freePar, double &cfac, vector<vector<size_t>> &iter, int &j, int &k) {
    double minSig = 1e-6;
    double eps=0.03;
    cfac *= (1+eps);
    if (cfac > 1e20) cfac = 1.0;

    // mixAllModels[0]: weight * par^2
    // mixAllModels[1]: weight
    // mixAllModels[2]: weight * par
    // mixAllModels[3]: new parameter
    // mixAllModels[4]: original parameter
    // mixAllModels[5]: best parameter
    // mixAllModels[6]: lower bound
    // mixAllModels[7]: upper bound
    // mixAllModels[8]: step
    j = random() * param.nLens;
    k = freePar[j][(int) (random() * freePar[j].size())];
    //cout << "parameter " << j << " " << k << endl;

    double stepSig = param.mixAllModels[8][j].paraList[k];
    double par0 = param.mixAllModels[4][j].paraList[k];

    param.mixAllModels[0][j].paraList[k] += cfac*par0*par0;
    param.mixAllModels[1][j].paraList[k] += cfac;
    param.mixAllModels[2][j].paraList[k] += cfac*par0;
    if (iter[j][k] > 3) {
        double sig = sqrt((1+1e-7)*param.mixAllModels[0][j].paraList[k]/param.mixAllModels[1][j].paraList[k]
                - pow(param.mixAllModels[2][j].paraList[k]/param.mixAllModels[1][j].paraList[k],2));
        //stepSig = 3.*2.38*sig/sqrt(freePar[j].size());
        stepSig = sig;
        if (std::isnan(stepSig)) {
            iter[j][k] = 0;
            param.mixAllModels[0][j].paraList[k] = 0.;
            param.mixAllModels[1][j].paraList[k] = 0.;
            param.mixAllModels[2][j].paraList[k] = 0.;
            stepSig = param.mixAllModels[8][j].paraList[k];
        }
    }
    if (stepSig < minSig) stepSig = minSig;
    double r = cgauss();
    param.mixAllModels[3][j].paraList[k] = par0 + r*stepSig;
    //std::cout << r << " " << stepSig << " " << std::endl;
    if (param.mixAllModels[3][j].paraList[k] < param.mixAllModels[6][j].paraList[k]
            || param.mixAllModels[3][j].paraList[k] > param.mixAllModels[7][j].paraList[k])
        param.mixAllModels[3][j].paraList[k] = param.mixAllModels[6][j].paraList[k] +
            random()*(param.mixAllModels[7][j].paraList[k]-param.mixAllModels[6][j].paraList[k]);

    iter[j][k]++;
    return stepSig;
}

void MC::stepPar(MultModelParam &param, vector<vector<size_t>> &freePar, double &cfac, int &iter) {
    double minSig = 1e-6;
    double eps=0.03;
    cfac *= (1+eps);
    if (cfac > 1e20) cfac = 1.0;

    for(int j=0; j<param.nLens; ++j) {
        for (auto k: freePar[j]) {
            double stepSig = param.mixAllModels[8][j].paraList[k];
            double par0 = param.mixAllModels[4][j].paraList[k];
            param.mixAllModels[0][j].paraList[k] += cfac*par0*par0;
            param.mixAllModels[1][j].paraList[k] += cfac;
            param.mixAllModels[2][j].paraList[k] += cfac*par0;
            if (iter > 3) {
                double sig = sqrt((1+1e-7)*param.mixAllModels[0][j].paraList[k]/param.mixAllModels[1][j].paraList[k]
                        - pow(param.mixAllModels[2][j].paraList[k]/param.mixAllModels[1][j].paraList[k],2));
                stepSig = sig/sqrt(freePar[j].size());
                if (std::isnan(stepSig)) {
                    iter = 0;
                    param.mixAllModels[0][j].paraList[k] = 0.;
                    param.mixAllModels[1][j].paraList[k] = 0.;
                    param.mixAllModels[2][j].paraList[k] = 0.;
                    stepSig = param.mixAllModels[8][j].paraList[k];
                }
            }
            if (stepSig < minSig) stepSig = minSig;
            double r = cgauss();
            param.mixAllModels[3][j].paraList[k] = par0 + r*stepSig;
            //std::cout << r << " " << stepSig << " " << std::endl;
            if (param.mixAllModels[3][j].paraList[k] < param.mixAllModels[6][j].paraList[k]
                    || param.mixAllModels[3][j].paraList[k] > param.mixAllModels[7][j].paraList[k])
                param.mixAllModels[3][j].paraList[k] = param.mixAllModels[6][j].paraList[k] +
                    random()*(param.mixAllModels[7][j].paraList[k]-param.mixAllModels[6][j].paraList[k]);
        }
    }
    iter++;
}

double MC::stepPar(vector<vec> &src, double cfac, size_t &iter) {
    double minSig = 1e-6;
    double eps=0.0;
    cfac *= (1+eps);
    if (cfac > 1e20) {
        cfac = 1.0;
        for (int k=0; k<src[0].size(); ++k) {
            src[0](k) = 0.;
            src[1](k) = 0.;
            src[2](k) = 0.;
        }
    }

    for (int k=0; k<src[0].size(); ++k) {
        double stepSig = src[8](k);
        double par0 = src[4](k);

        src[0](k) += cfac*par0*par0;
        src[1](k) += cfac;
        src[2](k) += cfac*par0;
        if (iter > 3) {
            double sig = sqrt(src[0](k)/src[1](k)
                    - pow(src[2](k)/src[1](k),2));
            stepSig = 3.*2.38*sig/sqrt(src[0].size());
            if (std::isnan(stepSig)) {
                iter=0;
                src[0](k) = 0.;
                src[1](k) = 0.;
                src[2](k) = 0.;
                stepSig = src[8](k);
            }
        }
        if (stepSig < minSig) stepSig = minSig;
        double r = cgauss();
        src[3][k] = par0 + r*stepSig;
        if (src[3](k) < src[6](k))
            src[3](k) = src[6](k);
        if (src[3](k) > src[7](k))
            src[3](k) = src[7](k);
    }
    return cfac;
}
