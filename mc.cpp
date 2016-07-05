/*
 * mc.cpp
 *
 *  Created on: Apr 23, 2016
 *      Author: En-Hsin Peng
 */

#include "mc.h"
#include <iostream>
#include <iomanip>

void MC::setRng(unsigned seed) {
    rng_engine.seed(seed);
    rng = std::bind(std::uniform_real_distribution<double>(0.,1.), std::ref(rng_engine));
    normal = std::bind(std::normal_distribution<double>(0.,1.), std::ref(rng_engine));
    //makeCgauss();
}

void MC::setParam(MultModelParam &param) {
    freePar.resize(param.nLens);
    iter.resize(param.nLens);
    nFreePar=0;
    for(int j=0; j<param.nLens; ++j) {
        iter[j].resize(param.mixAllModels[0][j].paraList.size(),0);
        for (size_t k=0; k<param.mixAllModels[0][j].paraList.size(); ++k) {
            if (param.mixAllModels[6][j].paraList[k] < param.mixAllModels[7][j].paraList[k]) {
                param.mixAllModels[4][j].paraList[k] = param.mixAllModels[6][j].paraList[k] +
                                random()*(param.mixAllModels[7][j].paraList[k]-param.mixAllModels[6][j].paraList[k]);
                freePar[j].push_back(k);
                nFreePar++;
            }
        }
    }
    cout << "total par " << nFreePar<< endl;
    iters=0;
    nAccept=0;
}
MC::MC(MultModelParam &param, unsigned seed) {
    setRng(seed);
    cfac=1.0;
    setParam(param);
}

MC::MC(Model* model, Conf* conf,
     const std::function<double(Model*)> &objective, size_t n, string outputFileName, size_t &iter) {
    setRng(conf->seed);
    setParam(model->param);
    setupGW(model->param, n);
    if (conf->GA) {
        setupGA();
        chkptFileName = "mcga_chkpt_"+to_string(conf->seed)+".txt";
    } else {
        chkptFileName = "mcgw_chkpt_"+to_string(conf->seed)+".txt";
    }
    if (conf->resume) {
        load(iter);
        output.open(outputFileName, std::ofstream::out | std::ofstream::app);
    } else {
        output.open(outputFileName);
    }
    this->objective = objective;
    if (conf->GA) setupGA();
}

void MC::setupGW(MultModelParam &param, size_t n) {
    nWalker = n;
    stepA = 2.;
    par.resize(n);
    par[0].resize(nFreePar);
    bound.resize(2);
    bound[0].resize(nFreePar);
    bound[1].resize(nFreePar);
    bestPar.resize(nFreePar,0);
    R0.resize(n,std::numeric_limits<double>::max());
    RMin = std::numeric_limits<double>::max();
    size_t c(0);
    double ss(5e-3);
    for (size_t j=0; j<freePar.size(); ++j) {
        for (auto k: freePar[j]) {
            bound[0][c] = param.mixAllModels[6][j].paraList[k];
            bound[1][c] = param.mixAllModels[7][j].paraList[k];
            par[0][c] = 0.5 * (bound[0][c] + bound[1][c]);
            par[0][c] += (random()-0.5) * 0.5 * (bound[1][c]-bound[0][c]);
            c++;
        }
    }
    for (size_t m=1; m<nWalker; ++m) {
        par[m].resize(nFreePar);
        for (size_t c=0; c<nFreePar; ++c) par[m][c] = par[m-1][c] + ss*(random()-0.5)*(bound[1][c]-bound[0][c]);
    }
}

void MC::setupGA() {
    pCrossOver = 0.85;
    pMutation = 0.05;
    index.resize(nWalker);
    for (size_t m=0; m<nWalker; ++m) {
        index[m] = m;
        for (size_t c=0; c<nFreePar; ++c) par[m][c] = bound[0][c] + random()*(bound[1][c]-bound[0][c]);
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

void MC::stretchMove(Model *model, size_t kk) {
    size_t jj;
    do {
        jj = random()*nWalker;
    } while (jj == kk);

    size_t c(0);
    double zz = pow((stepA - 1)*random() + 1,2)/stepA;
    for(int j=0; j<model->param.nLens; ++j) {
        for (auto k: freePar[j]) {
            model->param.mixAllModels[3][j].paraList[k] = par[jj][c] + zz*(par[kk][c]-par[jj][c]);
            if (model->param.mixAllModels[3][j].paraList[k]<bound[0][c] || model->param.mixAllModels[3][j].paraList[k]>bound[1][c])
                return;
            c++;
        }
    }
    iters++;
    double R = objective(model);
    if (R < 0) return;
    double weight(0.5);
    if (R < R0[kk] || log(random()) <= (nFreePar-1)*log(zz) - (R-R0[kk])*weight) {
        R0[kk] = R;
        size_t c(0);
        for(int j=0; j<model->param.nLens; ++j) {
            for (auto k: freePar[j]) {
                par[kk][c] = model->param.mixAllModels[3][j].paraList[k];
                if (R < RMin) bestPar[c] = par[kk][c];
                c++;
            }
        }
        nAccept++;
    }
    if (R < RMin) RMin = R;
}

void MC::load(size_t &loop) {
    ifstream input(chkptFileName.c_str());
    string line, token;
    if (input) {
        for (size_t j=0; j<nWalker; ++j) {
            getline(input, line);
            istringstream ss(line);
            for (size_t m=0; m<nFreePar; ++m) {
                getline(ss, token, ' ');
                par[j][m] = stod(token);
            }
            getline(ss, token, ' ');
            R0[j] = stod(token);
        }
        getline(input, line);
        istringstream ss(line);
        for (size_t m=0; m<nFreePar; ++m) {
            getline(ss, token, ' ');
            bestPar[m] = stod(token);
        }
        getline(ss, token, ' ');
        RMin = stod(token);
        getline(input, line);
        loop = stoi(line);
    }
}



void MC::checkPoint(size_t loop) {
    ofstream chkpt;
    chkpt.open(chkptFileName);
    chkpt << std::scientific << std::setprecision(4);
    for (size_t j=0; j<nWalker; ++j) {
        for (size_t m=0; m<nFreePar; ++m) {
            chkpt << par[j][m] << " ";
        }
        chkpt << R0[j] << endl;
    }
    for (size_t m=0; m<nFreePar; ++m) {
        chkpt << bestPar[m] << " ";
    }
    chkpt << RMin << endl;
    chkpt << loop << endl;
    chkpt.close();
}

void MC::writeOutput(size_t loop) {
    output << std::scientific << std::setprecision(4);
    for (size_t m=0; m<nFreePar; ++m) {
        output << bestPar[m] << " ";
    }
    output << RMin << " ";
    if (iters>0) output << std::fixed << (double) nAccept/iters << " ";
    output << std::setw(10) << loop << endl;
    //output.close();  //let destructor close it
}

void MC::writeOutput(size_t loop, int thin) {
    double rate = (double) nAccept/iters;
    for (size_t k=0; k<nWalker/thin; ++k) {
        size_t j = random()*nWalker;
        output << std::scientific << std::setprecision(4);
        for (size_t m=0; m<nFreePar; ++m) {
            output << par[j][m] << " ";
        }
        output << R0[j] << " ";
        output << std::fixed << rate << " ";
        output << std::setw(10) << loop << endl;
    }
    nAccept = 0;
    iters = 0;
}

void MC::copyParam(MultModelParam &param) {
    size_t c(0);
    for(int j=0; j<param.nLens; ++j) {
        for (auto k: freePar[j]) {
            param.mixAllModels[5][j].paraList[k] = bestPar[c];
            c++;
        }
    }
}


void MC::evaluate(Model *model, size_t m) {
    size_t c(0);
    for(int j=0; j<model->param.nLens; ++j) {
        for (auto k: freePar[j]) {
            model->param.mixAllModels[3][j].paraList[k] = par[m][c];
            c++;
        }
    }
    R0[m] = objective(model);
}

void MC::startGA() {
    //selection
    vector<vector<double>> parPrev(par);
    for (size_t m=0; m<nWalker; ++m) {
        size_t jj = pow(random(),1.7) * nWalker;
        for (size_t k=0; k<nFreePar; ++k) par[m][k] = parPrev[index[jj]][k];
    }

    //crossover
    int first(0);
    size_t i;
    for (size_t m=0; m<nWalker; ++m) {
        double p = random();
        if (p < pCrossOver) {
            ++first;
            if (first & 1) {  //odd
                i = m;
            } else {
                for (size_t k=0; k<random()*nFreePar; ++k) {
                    size_t kk = random()*nFreePar;
                    double t = par[i][kk];
                    par[i][kk] = par[m][kk];
                    par[m][kk] = t;
                }
            }
        }
    }

    //mutation
    for (size_t t=0; t<pMutation*nFreePar*nWalker; ++t) {
        int p = random()*nFreePar*nWalker;
        size_t j = p / nFreePar;
        size_t k = p % nFreePar;
        par[j][k] = bound[0][k] + random()*(bound[1][k] - bound[0][k]);
    }
}

void MC::elitism() {
    index = sort_indexes(R0);
    if (R0[index[0]] < RMin) {
        RMin = R0[index[0]];
        for (size_t k=0; k<nFreePar; ++k) bestPar[k] = par[index[0]][k];
    } else {
        for (size_t k=0; k<nFreePar; ++k) par[index[0]][k] = bestPar[k];
        R0[index[0]] = RMin;
    }
}

