#include "statistics.h"

void test_all(int max_ind, std::string folder) {
	std::string range_str = "def_range.alv";
	std::string E_str = "E_wcrit_distro.alv";
	std::string yield_str = "yield_wcrit_distro.alv";

	for (int i = 0; i <= max_ind; i++) {

		//считаем данные
		std::string range_str_i = folder + "/" + std::to_string(i) + range_str;
		std::string E_str_i = folder + "/" + std::to_string(i) + E_str;
		std::string yield_str_i = folder + "/" + std::to_string(i) + yield_str;

		indexData E_data, yield_data;
		double range;

		int num_E;
		std::fstream  file_E(E_str_i);
		while (file_E >> num_E)
			E_data.push_back(num_E);

		int num_yield;
		std::fstream  file_yield(yield_str_i);
		while (file_yield >> num_yield)
			yield_data.push_back(num_yield);

		std::fstream  file_range(range_str_i);
		file_range >> range;

		//обрабатываем
		intPointVec grouped_data_E, grouped_data_yield;
		group_data(grouped_data_E, E_data);
		group_data(grouped_data_yield, yield_data);

		double E_mean, yield_mean;
		E_mean = calc_mean(grouped_data_E, E_data);
		yield_mean = calc_mean(grouped_data_yield, yield_data);

		int3PointVec observed_expected_E, observed_expected_yield;
		calc_expected(observed_expected_E, E_mean, grouped_data_E, E_data);
		calc_expected(observed_expected_yield, yield_mean, grouped_data_yield, yield_data);

		double chi_square_E, chi_square_yield;
		chi_square_E =  calc_chi_square(observed_expected_E);
		chi_square_yield = calc_chi_square(observed_expected_yield);

		std::cout << "------------------------------------" << std::endl;
		std::cout << "test " << i << std::endl;
		std::cout << "range = " << range << ";" << std::endl;
		std::cout << "degrees = " << grouped_data_yield.size()-1 << std::endl;
		std::cout << "E_mean = " << E_mean << "; E_chi_square = " << chi_square_E << ";" << std::endl;
		std::cout << "yield_mean = " << yield_mean << "; yield_chi_square = " << chi_square_yield << ";" << std::endl;
		std::cout << "------------------------------------" << std::endl;
	}
}

void group_data(intPointVec& grouped_data, const indexData& data) {

	grouped_data.push_back(intPoint(data[0], 1));

	for (int i = 1; i < data.size(); i++) {
		bool found = false;
		for (int j = 0; j < grouped_data.size(); j++) 
			if (grouped_data[j].val == data[i]) {
				grouped_data[j].occur += 1;
				found = true;
			}

		if(!found)
			grouped_data.push_back(intPoint(data[i], 1));
	}
}

double calc_mean(const intPointVec& grouped_data, const indexData& data) {
	double val = 0;
	for(int i = 0; i < grouped_data.size(); i++)
		val += (grouped_data[i].occur * grouped_data[i].val);
	val /= data.size();

	return val;
}

double geometric_density(double mean, int n) {
	double p = 1.0 / mean;
	return pow((1.0 - p), n - 1) * p;
}

void calc_expected(int3PointVec& observed_expected, double mean, const intPointVec& grouped_data, const indexData& data) {
	for (int i = 0; i < grouped_data.size(); i++) {
		double prob = geometric_density(mean, grouped_data[i].val) * data.size();
		observed_expected.push_back(int3Vec(grouped_data[i].val, grouped_data[i].occur, prob));
	}
}

double calc_chi_square(const int3PointVec& observed_expected) {
	double val = 0;
	for (int i = 0; i < observed_expected.size(); i++)
		val += (pow(observed_expected[i].observed - observed_expected[i].expected, 2) 
		/ observed_expected[i].expected);

	return val;
}