#include "helper.h"

vector<double> GenerateVector(int n)
{
	vector<double> vector(n);
	for (int i = 0; i < n; i++)
	{
		vector[i] = ((double)rand()) / RAND_MAX * 2 - 1; //rand()%1+1;
	}
	return vector;
}
//vector[i] = ((double)rand()) / RAND_MAX * 2 - 1;