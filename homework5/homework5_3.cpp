#include <iostream>
#include <cmath>
#include <random>
#include <iomanip>
#include <vector>
#include <stack>
#include <queue>

using namespace std;

#include "../Matrix.hpp"
#include "../Interpolation.hpp"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
	Interpolation Runge_Interpolation(11);
	Runge_Interpolation.Shooting_method(101);
	// Runge_Interpolation.Print_Interpolation_y();

	Runge_Interpolation.B_spline(101);
	Runge_Interpolation.Print_Interpolation_y();


	return 0;
}