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
    MC(MultModelParam &param, unsigned seed, size_t n);
    double gauss(double sigma, double mu, double x);
    void makeCgauss();
    double cgauss();
    double random();
    double stepPar(MultModelParam &param, int &j, int &k);
    void stepPar(MultModelParam &param);
    double stepPar(vector<vec> &src, double cfac, size_t &iter);
    void checkPoint(string outputFileName, MultModelParam &param, bool moveAll, double L, double LMax);
    void checkPoint(string outputFileName, MultModelParam &param, vector<double> &R0, double RMin, size_t loop);
    void load(string fileName, MultModelParam &param, bool moveAll, double &L, double &LMax);
    void load(string fileName, MultModelParam &param, vector<double> &R0, double &RMin, size_t &loop);
    void printIterNum(bool moveAll);
    double strechMove(MultModelParam &param, size_t kk);
    void updateMove(MultModelParam &param, size_t kk);
    void writeOutput(ofstream &output, MultModelParam &param, double RMin, size_t loop, double rate);
    void writeOutput(ofstream &output, vector<double> &R0, size_t loop, int thin, double rate);

private:
    std::vector<double> cgArr;
    std::mt19937 rng_engine;
    //std::default_random_engine rng_engine;
    std::function<double()> rng;
    std::function<double()> normal;
    double cfac;
    vector<vector<size_t>> freePar, iter;
    size_t iters;

    void setRng(unsigned seed);
    void setParam(MultModelParam &param);

    //GW
    size_t nFreePar;
    size_t nWalker;
    double stepA;
    vector<vector<double>> par, bound;
    void setupGW(MultModelParam &param, size_t n);
};

#endif /* MC_H_ */
