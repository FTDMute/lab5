#include "MyMath.h"
using namespace std;


RibbonMatrix MultiplicationRibbonMatrixOnTransponate(RibbonMatrix& mat)
{
	RibbonMatrix newMat = RibbonMatrix(mat.n, 2 * mat.l, 1);
	int n = newMat.n;
	RibbonMatrix mat2 = TransponateRibbonMatrix(mat);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			double koef = 0;
			for (int k = 0; k < n; k++)
			{
				koef += mat2.GetIJElemnet(i, k) * mat.GetIJElemnet(k, j);
			}
			newMat.SetIJElemnet(i, j, koef);
		}
	}


	return newMat;

}

RibbonMatrix TransponateRibbonMatrix(RibbonMatrix& mat)
{
	RibbonMatrix newMat = RibbonMatrix(mat.n, mat.l, 1);

	for (int i = 0; i < mat.n; i++)
	{
		for (int j = 0; j < mat.n; j++)
		{
			newMat.SetIJElemnet(j, i, mat.GetIJElemnet(i, j));
		}
	}

	return newMat;


}

vector<double> MultiplicateRibbonOnVector(RibbonMatrix& mat1, vector<double>& vec) // матрица слева, вектор справа
{
	vector<double> newVec(vec.size());
	int n = vec.size();
	for (int i = 0; i < n; i++)
	{
		double koef = 0;
		for (int k = 0; k < n; k++)
		{
			newVec[i] += mat1.GetIJElemnet(i, k) * vec[k];
		}
	}

	return newVec;
}

RibbonMatrix MultiplicatedMatrix(RibbonMatrix mat1, double ch) {
	RibbonMatrix MultiplicatedMatrix(mat1.n, mat1.l, 1);
	for (int i = 0; i < mat1.n; i++) {
		for (int j = 0; j < mat1.n; i++) {
			MultiplicatedMatrix.SetIJElemnet(i, j, mat1.GetIJElemnet(i, j) * ch);
		}
	}
	return MultiplicatedMatrix;
}

double ScalarMultiplication(vector<double> vec1, vector<double> vec2)
{
	double sum = 0;
	for (int i = 0; i < vec1.size(); i++)
		sum += vec1[i] * vec2[i];
	return sum;
}

vector<double> Minus(vector<double> vec, vector<double> vec1)
{
	vector<double> rez; rez = vec;
	for (int i = 0; i < vec.size(); i++)
		rez[i] = vec[i] - vec1[i];
	return rez;
}

vector<double> Plus(vector<double> vec, vector<double> vec1)//сумма векторов
{
	vector<double> rez; rez = vec;
	for (int i = 0; i < vec.size(); i++)
		rez[i] = vec[i] + vec1[i];
	return rez;
}

vector<double> MultiplicationNumberOnVetor(double ch, vector<double> vec)//Умножение число на вектор
{
	for (int i = 0; i < vec.size(); i++)
	{
		vec[i] *= ch;
	}
	return vec;
}

double Norm(vector<double> vec)//Норма вектора
{
	double sum = 0;
	for (int i = 0; i < vec.size(); i++)
		sum += vec[i] * vec[i];
	return sqrt(sum);
}

double errorCalc( vector<double>& trueAnsw, vector<double>& currentAnsw, double error)
{
	for (int i = 0; i < trueAnsw.size(); i++)
	{
		error = max(error, abs(trueAnsw[i] - currentAnsw[i]));
	}
		//std::cout << error << endl;
	return error;
}

void Ykobi(vector<double>& x, RibbonMatrix& A, vector<double>& b)
{
	int n = b.size();
	double error = 200.;
	//cout << error << endl;
	vector <double> xCurrent(n);
	vector <double> xLast(b);
	int k = 0;

	while (error > 1e-5) {
		error = 0.;
		for (int i = 0; i < n; i++) {
			double sum = 0;
			for (int j = 0; j < n; j++)
			{

				if (i != j)
				{
					sum += A.GetIJElemnet(i, j) * xLast[j];
				}
			}
			xCurrent[i] = (1. / A.GetIJElemnet(i, i)) * (b[i] - sum);
		}

		for (int i = 0; i < n; i++)
		{
			xLast[i] = xCurrent[i];
		}

		error = errorCalc( x, xCurrent,error);
		//cout << error << endl;
		k++;
	}
	std::cout << "for Yakobi needed error on " << k << " iteration" << endl;
}

vector <double> YkobiForPcgm(vector<double>& x, RibbonMatrix& A, vector<double>& b)
{
	int n = b.size();
	double error = 200.;
	vector <double> xCurrent(n);
	vector <double> xLast(b);

	while (error > 1e-5) {
		error = 0.;
		for (int i = 0; i < n; i++) {
			double sum = 0;
			for (int j = 0; j < n; j++)
			{

				if (i != j)
				{
					sum += A.GetIJElemnet(i, j) * xLast[j];
				}
			}
			xCurrent[i] = (1. / A.GetIJElemnet(i, i)) * (b[i] - sum);
		}

		for (int i = 0; i < n; i++)
		{
			xLast[i] = xCurrent[i];
		}

		error = errorCalc(x, xCurrent, error);
	}
	return xCurrent;
}

void SOR(vector<double>& x, RibbonMatrix& A, vector<double>& b)
{
	int n = b.size();
	double error = 200.;
	//cout << error << endl;
	vector <double> xCurrent(n);
	vector <double> xLast(b);
	int k = 0;
	double x_k = 1;
	double w = 1;

	while (error > 1e-5) {
		error = 0.;
		for (int i = 0; i < n; i++) {
			double sum = 0.;
			double sum1 = 0.;
			if (i != 0) {
				for (int j = 0; j <= i-1; j++)
				{
					sum += A.GetIJElemnet(i, j) * xCurrent[j];
				}
			}
			for (int j = i+1; j < n; j++)
			{
					sum1 += A.GetIJElemnet(i, j) * xLast[j];
			}
			xCurrent[i] = (1. - w) * x_k + (w / A.GetIJElemnet(i, i)) * (b[i] - sum - sum1);
			error = max(abs(xCurrent[i]-xLast[i]),error);
			//cout << error << endl;
		}

		for (int i = 0; i < n; i++)
		{
			xLast[i] = xCurrent[i];
		}

		
		//cout << error << endl;
		k++;
	}
	std::cout << "for SOR needed error on " << k << " iteration" << endl;
}

void CGM(vector<double>& x, RibbonMatrix& A, vector<double>& b) {
	vector <double> x_0(b.size(), 1.);
	vector <double> A_x_0 = MultiplicateRibbonOnVector(A, x_0);
	vector <double> r = Minus(b, A_x_0);
	vector <double> r_0;
	vector <double> z = r;
	double a;
	double c;
	int k = 0;
	do {
		k++;
		a = ScalarMultiplication(r, r) / ScalarMultiplication(MultiplicateRibbonOnVector(A, z), z);
		x_0 = Plus(x_0, MultiplicationNumberOnVetor(a, z));
		r_0 = r;
		r = Minus(r, MultiplicationNumberOnVetor(a, MultiplicateRibbonOnVector(A, z)));
		c = ScalarMultiplication(r, r) / ScalarMultiplication(r_0, r_0);
		z = Plus(r, MultiplicationNumberOnVetor(c, z));
	} while (Norm(z) > 0.0001);
	std::cout << "for CGM needed error on " << k << " iteration" << endl;
}

void PCGM(vector<double>& x, RibbonMatrix& A, vector<double>& b) {
	vector <double> x_0=YkobiForPcgm(x,A,b);
	vector <double> A_x_0 = MultiplicateRibbonOnVector(A, x_0);
	vector <double> r = Minus(b, A_x_0);
	vector <double> r_0;
	vector <double> z = r;
	double a;
	double c;
	int k = 0;
	do {
		k++;
		a = ScalarMultiplication(r, r) / ScalarMultiplication(MultiplicateRibbonOnVector(A, z), z);
		x_0 = Plus(x_0, MultiplicationNumberOnVetor(a, z));
		r_0 = r;
		r = Minus(r, MultiplicationNumberOnVetor(a, MultiplicateRibbonOnVector(A, z)));
		c = ScalarMultiplication(r, r) / ScalarMultiplication(r_0, r_0);
		z = Plus(r, MultiplicationNumberOnVetor(c, z));
	} while (Norm(z) > 0.0001);
	std::cout << "for PCGM needed error on " << k << " iteration" << endl;
}