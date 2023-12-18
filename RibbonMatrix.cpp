#include "RibbonMatrix.h"

RibbonMatrix::RibbonMatrix(int n, int l, double q)
{
	this->n = n;
	this->l = l;
	GenerateAllRibbons();
	CalculateMainDiagonal(q);
}





double RibbonMatrix::GetIJElemnet(int i, int j)
{
	if (i == j)
	{
		return mainDiagonal[i];
	}
	if (abs(j - i) - 1 < l)
	{
		if (j > i)
		{
			return upperRibbons[j - i - 1][i];
		}
		if (j < i)
		{
			return lowerRibbons[i - j - 1][j];
		}
	}
	return 0;
}

void RibbonMatrix::SetIJElemnet(int i, int j, double value)
{
	if (i == j)
		mainDiagonal[i] = value;
	else if (abs(j - i) - 1 < l)
	{
		if (j > i)
		{
			upperRibbons[j - i - 1][i] = value;
		}
		if (j < i)
		{
			lowerRibbons[i - j - 1][j] = value;
		}
	}
	else
	{
		//cout << " chto to ne tak";
	}


}


void RibbonMatrix::GenerateAllRibbons()
{
	upperRibbons = vector < vector<double>>(l);
	lowerRibbons = vector < vector<double>>(l);
	for (int i = 0; i < l; i++)
	{
		
		upperRibbons[i] = GenerateVector(n - i);
		lowerRibbons[i] = GenerateVector(n - i);

	}
}


void RibbonMatrix::CalculateMainDiagonal(double q)
{
	mainDiagonal = vector<double>(n, 0);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < l; j++)
		{
			if (i < upperRibbons[j].size())
			{
				mainDiagonal[i] += fabs(upperRibbons[j][i]);
			}
			if (i - 1 - j < lowerRibbons[j].size() && i - 1 - j >=0)
			{
				mainDiagonal[i] += fabs(lowerRibbons[j][i - 1 - j]);

			}
		}
		mainDiagonal[i] *= q;

	}
}
