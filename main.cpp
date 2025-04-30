//Center difference Newton method based on variable step size without backward difference
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream> 
#include <string>
using namespace std;

double Gamma = 1.4;//Ratio of Specific Heat
double G_add = Gamma + 1;
double G_subtract = Gamma - 1;

//Laval nozzle formula
double f(double Ma, double A) {
	if (Ma == 0)
	{
		return 1.0;
	}
	double Ma2 = pow(Ma, 2);
	double temp = pow(2 / G_add * (1 + G_subtract / 2 * Ma2), G_add / G_subtract) / Ma2 - pow(A, 2);
	return temp;
}

int main() {
	//Define Parameters
	double L;//nozzle length
	double dx;//Discrete distance
	int N;//Number of discrete points
	N = 31;
	vector<double> x(N, 0), A(N, 0), Ma(N, 1);//Distance Function

	L = 3.0;
	dx = L / (N - 1);
	for (int i = 0; i < N; i++)
	{
		x[i] = dx * i;
		A[i] = 1 + 2.2 * pow(x[i] - 1.5, 2);
	}

	/*	Ma = {
		 0.098,0.11,0.124,0.14,0.16,0.185,0.214,0.249,0.293,0.347,
		 0.413,0.494,0.592,0.709,0.846,1,1.169,1.348,1.531,1.715,1.896,
		 2.071,2.24,2.402,2.557,2.706,2.848,2.983,3.114,3.239,3.359
		};Analytical solution with three decimal places retained*/

		//Solve point by point
	ofstream outFile("C:\\Desktop\\CDNM.txt");
	for (int i = 0; i < N; i++)
	{
		double temp = f(Ma[i], A[i]);
		double temp_b, temp_f;
		double dMa = 0.1;
		temp_b = f(Ma[i] - dMa, A[i]);
		//Iterative solution
		while (abs(temp) > 1e-13) {
			temp_f = f(Ma[i] + dMa, A[i]);
			dMa = -2 * temp * dMa / (temp_f - temp_b);
			Ma[i] = Ma[i] + dMa;//Ma_New-dMa=Ma
			temp_b = temp;//temp_b=f(Ma_New-dMa,A)=f(Ma,A)=temp
			temp = f(Ma[i], A[i]);
		}
		outFile << setw(8) << x[i] << setw(16) << Ma[i] << setw(18) << temp << "\n";
	}
	outFile.close();

	return 0;
}