//Center difference Newton method based on variable step size that does not require backward calculation of function values
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream> 
#include <string>
using namespace std;

const double L = 3;//nozzle length
const int N = 31;//Number of discrete points
const double dx = L / (N - 1);//Discrete distance
const double residual = 1e-13;//Iterative accuracy
const double stepL = 0.01;//Initial step size
const double IMa = 0;//Initial value

//Laval nozzle formula
double f(const int i, double Ma) {
	if (Ma == 0)
	{
		return std::numeric_limits<double>::max();
	}
	double x, A;//Distance Function
	x = dx * i;
	A = 1 + 2.2 * pow(x - 1.5, 2);
	double Gamma = 1.4;//Ratio of Specific Heat
	double G_add = Gamma + 1;
	double G_subtract = Gamma - 1;
	double temp = pow(2 / G_add * (1 + G_subtract / 2 * pow(Ma, 2)), G_add / G_subtract / 2) / Ma - A;
	return temp;
}

void ForwardDifference(vector<double> Ma) {
	int j = 0;//The number of times the iterative process calculates the function(f)
	cout << "********************************************" << endl;
	for (int i = 0; i < N; i++)
	{
		double temp = f(i, Ma[i]);
		double temp_f;
		double dMa = stepL;
		//Iterative solution
		while (abs(temp) > residual) {
			temp_f = f(i, Ma[i] + dMa);
			dMa = -temp * dMa / (temp_f - temp);
			if (abs(dMa) < residual)
			{
				dMa = residual;
			}
			if (abs(dMa) > 1 / residual)
			{
				dMa = stepL;
			}
			Ma[i] = Ma[i] + dMa;
			temp = f(i, Ma[i]);
			++j;
		}
		cout << setw(8) << dx * i << setw(16) << Ma[i] << setw(18) << temp << endl;
	}
	cout << setw(8) << "Method" << setw(16) << "Times of f" << endl;
	cout << setw(8) << "FD:" << setw(16) << 2 * j << endl;
	cout << "********************************************" << endl;
}

void CenterDifference(vector<double> Ma) {
	int j = 0;//The number of times the iterative process calculates the function(f)
	cout << "********************************************" << endl;
	for (int i = 0; i < N; i++)
	{
		double temp = f(i, Ma[i]);
		double temp_b, temp_f;
		double dMa = stepL;
		if (abs(temp) > residual)
		{
			temp_b = f(i, Ma[i] - dMa);
		}
		//Iterative solution
		while (abs(temp) > residual) {
			temp_f = f(i, Ma[i] + dMa);
			dMa = -temp * 2 * dMa / (temp_f - temp_b);
			if (abs(dMa) < residual)
			{
				dMa = residual;
			}
			if (abs(dMa) > 1 / residual)
			{
				dMa = stepL;
			}
			Ma[i] = Ma[i] + dMa;//Ma_New=Ma+dMa_New
			temp_b = temp;//temp_b_New=f(Ma_New-dMa_New,A)=f(Ma,A)=temp
			temp = f(i, Ma[i]);
			++j;
		}
		cout << setw(8) << dx * i << setw(16) << Ma[i] << setw(18) << temp << endl;
	}
	cout << setw(8) << "Method" << setw(16) << "Times of f" << endl;
	cout << setw(8) << "CD:" << setw(16) << 2 * j << endl;
	cout << "********************************************" << endl;
}

void BackwardDifference(vector<double> Ma) {
	int j = 0;//The number of times the iterative process calculates the function(f)
	cout << "********************************************" << endl;
	for (int i = 0; i < N; i++)
	{
		double temp = f(i, Ma[i]);
		double temp_b;
		double dMa = stepL;
		if (abs(temp) > residual)
		{
			temp_b = f(i, Ma[i] - dMa);
		}
		//Iterative solution
		while (abs(temp) > residual) {
			dMa = -temp * dMa / (temp - temp_b);
			if (abs(dMa) < residual)
			{
				dMa = residual;
			}
			if (abs(dMa) > 1 / residual)
			{
				dMa = stepL;
			}
			Ma[i] = Ma[i] + dMa;//Ma_New=Ma+dMa_New
			temp_b = temp;//temp_b_New=f(Ma_New-dMa_New,A)=f(Ma,A)=temp
			temp = f(i, Ma[i]);
			++j;
		}
		cout << setw(8) << dx * i << setw(16) << Ma[i] << setw(18) << temp << endl;
	}
	cout << setw(8) << "Method" << setw(16) << "Times of f" << endl;
	cout << setw(8) << "BD:" << setw(16) << j << endl;
	cout << "********************************************" << endl;
}

void CFD(vector<double> Ma) {
	int j = 0;//The number of times the iterative process calculates the function(f)
	int maxB = 0;
	int sumB = 0;
	cout << "********************************************" << endl;
	for (int i = 0; i < N; i++)
	{
		double temp = f(i, Ma[i]);
		double temp_b, temp_f;
		double dMa = stepL;
		if (abs(temp) > residual)
		{
			maxB = static_cast<int>(std::round(log(abs(temp)) - log(residual)));//|temp|/e^(maxB)=residual
			maxB = static_cast<int>(std::round(0.5 * pow(1 + 8 * maxB, 0.5) - 0.5));//|temp|/e^[(1+maxB)maxB/2]=residual		
			sumB += maxB;
			maxB += j;
			temp_b = f(i, Ma[i] - dMa);
		}
		//Iterative solution
		while (abs(temp) > residual) {
			temp_f = f(i, Ma[i] + dMa);
			if (j < maxB)
			{
				dMa = -temp * 2 * dMa / (temp_f - temp_b);
			}
			else
			{
				dMa = -temp * dMa / (temp_f - temp);
			}
			if (abs(dMa) < residual)
			{
				dMa = residual;
			}
			if (abs(dMa) > 1 / residual)
			{
				dMa = stepL;
			}
			Ma[i] = Ma[i] + dMa;//Ma_New=Ma+dMa_New
			temp_b = temp;//temp_b_New=f(Ma_New-dMa_New,A)=f(Ma,A)=temp
			temp = f(i, Ma[i]);
			++j;
		}
		cout << setw(8) << dx * i << setw(16) << Ma[i] << setw(18) << temp << endl;
	}
	cout << setw(8) << "Method" << setw(16) << "Times of f" << setw(18) << "Times of f(x-dx)" << endl;
	cout << setw(8) << "CFD:" << setw(16) << 2 * j << setw(18) << sumB << endl;
	cout << "********************************************" << endl;
}

void BFD(vector<double> Ma) {
	int j = 0;//The number of times the iterative process calculates the function(f)
	int k = 0;
	int maxB = 0;
	int sumB = 0;
	cout << "********************************************" << endl;
	for (int i = 0; i < N; i++)
	{
		double temp = f(i, Ma[i]);
		double temp_b;
		double dMa = stepL;
		if (abs(temp) > residual)
		{
			maxB = static_cast<int>(std::round(log(abs(temp)) - log(residual)));//|temp|/e^(maxB)=residual
			maxB = static_cast<int>(std::round(0.5 * pow(1 + 8 * maxB, 0.5) - 0.5));//|temp|/e^[(1+maxB)maxB/2]=residual

			sumB += maxB;
			maxB += j;
			temp_b = f(i, Ma[i] - dMa);
		}
		//Iterative solution
		while (abs(temp) > residual) {
			if (j < maxB)
			{
				dMa = -temp * dMa / (temp - temp_b);
			}
			else
			{
				double temp_f = f(i, Ma[i] + dMa);
				dMa = -temp * dMa / (temp_f - temp);
				++k;
			}
			if (abs(dMa) < residual)
			{
				dMa = residual;
			}
			if (abs(dMa) > 1 / residual)
			{
				dMa = stepL;
			}
			Ma[i] = Ma[i] + dMa;//Ma_New=Ma+dMa_New
			temp_b = temp;//temp_b_New=f(Ma_New-dMa_New,A)=f(Ma,A)=temp
			temp = f(i, Ma[i]);
			++j;
		}
		cout << setw(8) << dx * i << setw(16) << Ma[i] << setw(18) << temp << endl;
	}
	cout << setw(8) << "Method" << setw(16) << "Times of f" << setw(18) << "Times of f(x-dx)" << setw(18) << "Times of f(x+dx)" << endl;
	cout << setw(8) << "BFD:" << setw(16) << j + k << setw(18) << sumB << setw(18) << k << endl;
	cout << "********************************************" << endl;
}

int main() {
	/*
		Ma = {
		 0.098,0.11,0.124,0.14,0.16,0.185,0.214,0.249,0.293,0.347,
		 0.413,0.494,0.592,0.709,0.846,1,1.169,1.348,1.531,1.715,1.896,
		 2.071,2.24,2.402,2.557,2.706,2.848,2.983,3.114,3.239,3.359
		};//Analytical solution with three decimal places retained
	*/

	//Initialize Ma
	vector<double> Ma(N, IMa);
	//The following section was not included during performance testing
  /*for (int i = 0; i < N; i++)
	{
		if (i > (N - 1) / 2)
		{
			Ma[i] = 2.0;
		}
	}*/

	cout << setw(8) << "x" << setw(16) << "Ma" << setw(18) << "temp" << endl;
	//Iterative solution
	ForwardDifference(Ma);
	CenterDifference(Ma);
	BackwardDifference(Ma);
	CFD(Ma);
	BFD(Ma);

	return 0;
}