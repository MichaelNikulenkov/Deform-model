#pragma once
#pragma once
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "nlopt.hpp"

struct Point {
	Point(double x = 0.0, double y = 0.0) : x(x), y(y) {};
	double x;
	double y;
};

typedef std::vector<std::vector<Point>> dataPoint;
typedef std::vector<std::vector<double>> dataDouble;
typedef std::vector<int> indexData;
typedef std::vector<double> vec;
typedef std::vector<std::vector<int>> dataInt;
typedef std::vector<Point> vecPoint;

struct CopyUtility {
	dataPoint transformed_points;
	dataPoint original_points;
	vec yield_stress;
	vec yield_deform;
	dataInt odquist_scale_ind;
	dataDouble odquist_intervals;
	vec odquist;
	vec d;
	vec modulus;
	double average_modulus;
	int current_interval_num;
};

void deleteEmptyLines(const std::string &FilePath);

//считать из файла точки, перевести в напр€жени€ и относ. деформации
//и преобразовать их в локальные ситемы координат циклов
void read_data(dataPoint& original_points, dataPoint& transformed_points, std::string filename);

void read_data_no_scale(dataPoint& original_points, dataPoint& transformed_points, std::string filename);

//модуль упругой деформации дл€ каждой точки в виде отношени€ сумм и усредненный модуль
void calc_current_modulus(const dataPoint& transformed_points, dataDouble& current_modulus);

//предел пропорцинальности
void calc_proportional_stress(const dataPoint& transformed_points, const dataDouble& current_modulus, indexData& pr_stress_ind);

//усредненный модуль упругости и модуль упругости дл€ каждого полуцикла
void calc_modulus(const dataDouble& current_modulus, const indexData& pr_stress_ind, vec& modulus, double& average_modulus);

//параметр одквиста
void calc_odqist(const dataPoint& original_points, vec& odquist);

//предел текучести
void calc_yield_stress(const dataPoint& transformed_points, const vec& modulus, const indexData& pr_stress_ind, vec& yield_stress, vec& yield_deform);

//разбиение полуциклов по шкале одквиста
void set_odquist_scale(const vec& odquist, dataDouble& odquist_intervals, dataInt& odquist_scale_ind);

//интерпол€ционный многочлен лагранжа дл€ первого полуцикла
double f(double deform, const dataPoint& transformed_points, double init_yield_deform, double init_modulus);

//производна€ полинома
double fd(double deform, const dataPoint& transformed_points, double init_yield_deform, double init_modulus);

//напр€жение по модели
double model_stress(double deform, double a, double b, double d, double yield_deform_init,
double yield_stress_init, double init_modulus, const CopyUtility& s);

double model_elastic_stress(double deform, double a, double b, double d, double yield_deform_init,
	double yield_stress_init, double init_modulus, const CopyUtility& s);

double model_stress_alpha(double deform, double a, double b, double d, double yield_deform_init,
	double yield_stress_init, double init_modulus, const CopyUtility& s, double alpha);

//функционал дл€ определени€ a,b. Cм. документацию NLopt. n = 2
double J(unsigned n, const double *x, double *grad, void* interval_num);

//расчет параметров a,b,d
void calc_d(const vec& modules, vec& d);
void calc_ab(double average_modulus, vec& a, vec& b, CopyUtility& s);

//интерпол€ционный полином дл€ параметров a, b, w
double interpolate_ab(double x, const vec& param, const vec& odquist_medium);
double interpolate_d(double x, const vec& d_vec, const vec& odquist, int halcycle_num);
//double interpolate_w(const vec& w, const dataInt& odquist_scale_ind, int halfcycle_ind);

//дл€ статистики
void calc_model_odqist_given_init_modulus(vec& model_odquist, double E_init, const CopyUtility& s, const vec& a, const vec& b, const vec& alpha, const vec& odquist_medium);
void calc_model_odqist_given_init_yield_stress(vec& model_odquist, double yield_deform_init, double yield_str_init, const CopyUtility& s, const vec& a, const vec& b, const vec& alpha, const vec& odquist_medium);
double calc_model_prop_and_yield_stress(double E_init, const CopyUtility& s);

void clear_file(std::string filename);
void write_to_file(const vec& x, const vec& y, std::string filename);
void write_to_file(int ind, const vec& x, const vec& y, std::string filename);
void write_to_file(const vec& x, std::string filename);
void write_to_file(const dataPoint& x, std::string filename);
void write_to_file(const std::vector<Point>& x, std::string filename);
void write_to_file(const dataDouble& x, std::string filename);
void write_to_file(const indexData& x, std::string filename);
void write_to_file_reverse(const std::vector<Point>& x, std::string filename);
void write_to_file_second_coord(const std::vector<Point>& x, std::string filename);

long double fact(int N);

double destruct_crit(double kappa, double gamma, double delta);
void calc_with_destruction_crit(indexData& destruct_ind_distro, dataDouble& model_odquist_destruct, const dataDouble& model_odquist);
void calc_with_destruction_crit_continuous(vec& destruct_ind_distro, dataDouble& model_odquist_destruct, const dataDouble& model_odquist);

double deformations_range(dataPoint& original_points);

double interpolate(double x, vecPoint points);

Point calc_criteria_point(std::string name);