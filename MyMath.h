#pragma once
#include <vector>
#include <iostream>
#include <ctime>
#include "RibbonMatrix.h"
using namespace std;



RibbonMatrix MultiplicationRibbonMatrixOnTransponate(RibbonMatrix& mat);
RibbonMatrix TransponateRibbonMatrix(RibbonMatrix& mat);
vector<double> MultiplicateRibbonOnVector(RibbonMatrix& mat1, vector<double>& vec);
vector<double> MultiplicateRibbonOnNumber(RibbonMatrix& mat1, double& ch, int n);
double errorCalc(vector<double>& trueAnsw, vector<double>& currentAnsw);
void Ykobi(vector<double>& x, RibbonMatrix& A, vector<double>& b);
void SOR(vector<double>& x, RibbonMatrix& A, vector<double>& b);
void CGM(vector<double>& x, RibbonMatrix& A, vector<double>& b);
void PCGM(vector<double>& x, RibbonMatrix& A, vector<double>& b);
//vector<double> SolutionOfSLAUWithLU(vector<vector<double>> a, vector<double> b);
//
//void Stroke(vector<vector<double>> a, vector<double>& x, vector<double>& y, int k);
//
//vector<double> ReverseStroke(vector<vector<double>> a, vector<double> y);
//
//
//vector<double> StaightStroke(vector<vector<double>> a, vector<double> y);
//void LU(vector<vector<double>> A, vector<vector<double>>& L, vector<vector<double>>& U, int n);

