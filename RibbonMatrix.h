#pragma once
#include <vector>
#include <math.h>
#include <iostream>
#include "helper.h"
using namespace std;
class RibbonMatrix
{
public:
	RibbonMatrix(int n, int l, double q); // генерирует случайную матрицу
	vector<vector<double>> upperRibbons;
	vector<vector<double>> lowerRibbons;
	vector<double> mainDiagonal;
	int l;
	int n;
	double GetIJElemnet(int i, int j);
	void SetIJElemnet(int i, int j, double value);

private:
	void GenerateAllRibbons(); 
	void CalculateMainDiagonal(double q);
};

