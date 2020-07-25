#pragma once
#include "func.h"
#include <vector>
#include <string>

struct intPoint {
	intPoint(int x = 0, int y = 0) : val(x), occur(y) {};
	int val;
	int occur;
};

struct int3Vec {
	int3Vec(int x = 0, int y = 0, double z = 0) : val(x), observed(y), expected(z) {};
	int val;
	int observed;
	double expected;
};

typedef std::vector<intPoint> intPointVec;
typedef std::vector<int3Vec> int3PointVec;

void group_data(intPointVec& grouped_data, const indexData& data);

double calc_mean(const intPointVec& grouped_data, const indexData& data);

double geometric_density(double mean, int n);

void calc_expected(int3PointVec& observed_expected, double mean, const intPointVec& grouped_data, const indexData& data);

double calc_chi_square(const int3PointVec& observed_expected);

void test_all(int max_ind, std::string folder);