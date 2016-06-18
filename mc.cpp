/*
 * mc.cpp
 *
 *  Created on: Apr 23, 2016
 *      Author: En-Hsin Peng
 */

#include "mc.h"
#include <iostream>
#include <iomanip>

MC::MC(MultModelParam &param, unsigned seed) {
    rng_engine.seed(seed);
    rng = std::bind(std::uniform_real_distribution<double>(0.,1.), std::ref(rng_engine));
    normal = std::bind(std::normal_distribution<double>(0.,1.), std::ref(rng_engine));
    //makeCgauss();
    cfac=1.0;
    iters=0;
    freePar.resize(param.nLens);
    iter.resize(param.nLens);
    for(int j=0; j<param.nLens; ++j) {
        iter[j].resize(param.mixAllModels[0][j].paraList.size(),0);
        for (size_t k=0; k<param.mixAllModels[0][j].paraList.size(); ++k) {
            if (param.mixAllModels[6][j].paraList[k] < param.mixAllModels[7][j].paraList[k]) {
                freePar[j].push_back(k);
            }
        }
    }
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


double MC::stepPar(MultModelParam &param, int &j, int &k) {
    double minSig = 1e-6;
    double eps=0.01;
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

void MC::stepPar(MultModelParam &param) {
    double minSig = 1e-6;
    double eps=0.01;
    cfac *= (1+eps);
    if (cfac > 1e20) cfac = 1.0;

    for(int j=0; j<param.nLens; ++j) {
        for (auto k: freePar[j]) {
            double stepSig = param.mixAllModels[8][j].paraList[k];
            double par0 = param.mixAllModels[4][j].paraList[k];
            param.mixAllModels[0][j].paraList[k] += cfac*par0*par0;
            param.mixAllModels[1][j].paraList[k] += cfac;
            param.mixAllModels[2][j].paraList[k] += cfac*par0;
            if (iters > 3) {
                double sig = sqrt((1+1e-7)*param.mixAllModels[0][j].paraList[k]/param.mixAllModels[1][j].paraList[k]
                        - pow(param.mixAllModels[2][j].paraList[k]/param.mixAllModels[1][j].paraList[k],2));
                stepSig = sig/sqrt(freePar[j].size());
                if (std::isnan(stepSig)) {
                    iters = 0;
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
    iters++;
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

void MC::checkPoint(string outputFileName, MultModelParam &param, bool moveAll, double L, double LMax) {
    ofstream chkpt;
    chkpt.open(outputFileName);
    for (int j=0; j<param.nLens; ++j) {
        for (auto k: freePar[j]) {
            for (size_t l=0; l<6; ++l) chkpt << std::scientific << std::setprecision(7) << param.mixAllModels[l][j].paraList[k] << " ";
            double par0 = param.mixAllModels[4][j].paraList[k];
            double s0 = param.mixAllModels[0][j].paraList[k] + cfac*par0*par0;
            double s1 = param.mixAllModels[1][j].paraList[k] + cfac;
            double s2 = param.mixAllModels[2][j].paraList[k] + cfac*par0;
            double sig = sqrt((1+1e-7)*s0/s1
                        - pow(s2/s1,2));
            chkpt << sig << " ";
            chkpt << endl;
        }
    }
    chkpt << cfac << endl;
    if (moveAll) {
        chkpt << iters << endl;
    } else {
        for (int j=0; j<param.nLens; ++j) {
            for (auto k: freePar[j]) {
                chkpt << iter[j][k] << " ";
            }
            chkpt << endl;
        }
    }
    chkpt << L << endl;
    chkpt << LMax << endl;
    chkpt.close();
}

void MC::load(string fileName, MultModelParam &param, bool moveAll, double &L, double &LMax) {
    ifstream input(fileName.c_str());
    string line, token;
    if (input) {
        for (int j=0; j<param.nLens; ++j) {
            for (auto k: freePar[j]) {
                getline(input, line);
                istringstream ss(line);
                for (size_t l=0; l<6; ++l) {
                    getline(ss, token, ' ');
                    param.mixAllModels[l][j].paraList[k] = stod(token);
                }
            }
        }
        getline(input, line);
        cfac = stod(line);
        if (moveAll) {
            getline(input, line);
            iters = stod(line);
        } else {
            for (int j=0; j<param.nLens; ++j) {
                if (freePar[j].size() > 0) {
                    getline(input, line);
                    istringstream ss(line);
                    for (auto k: freePar[j]) {
                        getline(ss, token, ' ');
                        iter[j][k] = stod(token);
                    }
                }
            }
        }
        getline(input, line);
        L = stod(line);
        LMax = stod(line);
        checkPoint("copy_"+fileName, param, moveAll, L, LMax);
    }
}

void MC::printIterNum(bool moveAll) {
    cout << "Iteration #: " << endl;
    if (moveAll) {
        cout << iters << endl;
    } else {
       for(size_t j=0; j<freePar.size(); ++j) {
           cout << "Lens " << j << endl;
           for (auto k: freePar[j]) {
               cout << iter[j][k] << " ";
           }
           cout << endl;
       }
    }
}
