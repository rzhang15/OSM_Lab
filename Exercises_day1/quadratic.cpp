#include <iostream>
#include <complex>
int main()
{
	using namespace std;
	cout << "Please enter the values of a, b, c: " <<endl;
	complex<float> a;
	complex<float> b;
	complex<float> c;
	cout <<"a = ";
	cin >> a; 
	cout <<"b = ";
	cin >> b;
	cout <<"c = ";
	cin >> c;
	complex<float> root1 = (-b+sqrt(b*b-float(4)*a*c))/(float(2)*a);
	complex<float> root2 = (-b-sqrt(b*b-float(4)*a*c))/(float(2)*a);
	cout << "Root 1 = " << root1 << endl;
	cout << "Root 2 = " << root2 << endl;
}

