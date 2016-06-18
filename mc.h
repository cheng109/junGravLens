/*
 * mc.h
 *
 *  Created on: Apr, 2016
 *      Author: En-Hsin Peng
 */

#ifndef MC_H_
#define MC_H_

#include <vector>
#include <random>
#include <functional>
#include <cmath>
#include "Model.h"

class MC {
public:
    MC(MultModelParam &param, unsigned seed);
    double gauss(double sigma, double mu, double x);
    void makeCgauss();
    double cgauss();
    double random();
    double stepPar(MultModelParam &param, int &j, int &k);
    void stepPar(MultModelParam &param);
    double stepPar(vector<vec> &src, double cfac, size_t &iter);
    void checkPoint(string outputFileName, MultModelParam &param, bool moveAll, double L, double LMax);
    void load(string fileName, MultModelParam &param, bool moveAll, double &L, double &LMax);
    void printIterNum(bool moveAll);

private:
    std::vector<double> cgArr;
    std::mt19937 rng_engine;
    //std::default_random_engine rng_engine;
    std::function<double()> rng;
    std::function<double()> normal;
    double cfac;
    vector<vector<size_t>> freePar, iter;
    size_t iters;
};

#endif /* MC_H_ */
