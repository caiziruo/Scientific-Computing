#ifndef Interpolation_hpp
#define Interpolation_hpp

#include <iostream>
#include <cmath>
#include <random>
#include <iomanip>
#include <vector>
#include <stack>
#include <queue>

using namespace std;

#include "Matrix.hpp"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Interpolation
{
public:
	int dimension_data;
	Matrix * data_t;
	Matrix * data_y;

	Interpolation(int n) {
		dimension_data = n;

		data_t = new Matrix(n, 1);
		data_y = new Matrix(n, 1);
	}

	void Newton_Interpolation();

	~Interpolation() {
		delete data_t;
		delete data_y;
	}
	
};

#endif