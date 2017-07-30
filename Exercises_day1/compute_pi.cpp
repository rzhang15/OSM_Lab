#include <iostream>
#include <math.h>

static long num_steps = 100000; 
double step;

int main()
{
	using namespace std;
	int i;
	double x, pi, sum = 0.0;
	step = 1.0/double(num_steps);
	for (i=0; i<num_steps; i++)
	{
		x = (i+0.5)*step;
		sum = sum + 4.0/(1.0+x*x);
	}
	pi = step*sum;
	cout << "The approximation of pi is: " << pi << endl;
	return 0;
}
