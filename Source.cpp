#include <iostream>
#include "helper.h"
#include "RibbonMatrix.h"
#include "MyMath.h"
#include <ctime>
using namespace std;

void GetAWithStarAndB(RibbonMatrix& A, vector<double>& b, vector<double>& x)
{
	b = MultiplicateRibbonOnVector(A, x);
	RibbonMatrix At = TransponateRibbonMatrix(A);
	RibbonMatrix AwithStar = MultiplicationRibbonMatrixOnTransponate( A);
	vector<double> bwithStar = MultiplicateRibbonOnVector(At, b);
	A = AwithStar;
	b = bwithStar;

}

void Task2(vector<double> x, RibbonMatrix A1, RibbonMatrix A2, RibbonMatrix A3, vector<double>b1, vector<double>b2, vector<double>b3) 
{
	cout << "Yakobi goes brrrrr" << endl;
	cout << "With q = 1.1 ";
	Ykobi(x, A1, b1);
	cout << "With q = 2 ";
	Ykobi(x,A2, b2);
	cout << "With q = 10 ";
	Ykobi(x, A3, b3);
	cout << endl;
}

void Task3(vector<double> x, RibbonMatrix A1, RibbonMatrix A2, RibbonMatrix A3, vector<double>b1, vector<double>b2, vector<double>b3)
{
	cout << "SOR goes brrrrr" << endl;
	cout << "With q = 1.1 ";
	SOR(x, A1, b1);
	cout << "With q = 2 ";
	SOR(x, A2, b2);
	cout << "With q = 10 ";
	SOR(x, A3, b3);
	cout << endl;
}

void Task4(vector<double> x, RibbonMatrix A1, RibbonMatrix A2, RibbonMatrix A3, vector<double>b1, vector<double>b2, vector<double>b3)
{
	cout << "CGM goes brrrrr" << endl;
	cout << "With q = 1.1 ";
	CGM(x, A1, b1);
	cout << "With q = 2 ";;
	CGM(x, A2, b2);
	cout << "With q = 10 ";
	CGM(x, A3, b3);
	cout << endl;
	cout << "PCGM goes brrrrr" << endl;
	cout << "With q = 1.1 " ;
	PCGM(x, A1, b1);
	cout << "With q = 2 " ;
	PCGM(x, A2, b2);
	cout << "With q = 10 " ;
	PCGM(x, A3, b3);
}


int main()
{
	srand(time(NULL));
	int n = 900, l = 25;
	int q1 = 1.1, q2 = 2, q3 = 10;
	int task;
	vector<double> x = GenerateVector(n);

	RibbonMatrix A1(n, l, q1);
	RibbonMatrix A2(n, l, q2);
	RibbonMatrix A3(n,l, q3);
	vector<double>b1;
	vector<double>b2;
	vector<double>b3;

	GetAWithStarAndB(A1, b1, x);
	GetAWithStarAndB(A2, b2, x);
	GetAWithStarAndB(A3, b3, x);

	cout << "enter task" << endl;
	cin >> task;
	//A1,A2,A3 и b1,b2,b3 нужные нам матрицы и они симметричные, положительно определенные, ленточные  и вектора соответсвенно
	switch (task) {
	case 1: {
		Task2(x, A1, A2, A3, b1, b2, b3);
		break;
	}
	case 2: {
		Task3(x, A1, A2, A3, b1, b2, b3);
		break;
	}
	case 3: {
		Task4(x, A1, A2, A3, b1, b2, b3);
	}
	}
	return 0;
}