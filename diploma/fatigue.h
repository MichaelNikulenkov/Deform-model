#pragma once
#include "halfcycle.h"
#include "nlopt.hpp"
//#include "nlopt.h"
#include <vector>
#include <memory>
#include "math.h"

//функционал для определения a,b. Cм. документацию NLopt, n = 2
double J(unsigned n, const double *x, double *grad, void* interval_num);

struct interval {
	//логарифмическая шкала до максимального значения накопленной пластической деформации
    std::vector<int> odqist_scale;
	//параметры a,b 
	std::vector<double> params{ 0.01, 0.01 };

	double odqist_init;
	double odqist_final;
};

//тк в nlopt нельзя минимизировать функции-члены-объекта, придется скопировать данные в отдельный объект
struct nlopt_copy_utility {
	int current_interval_num = 0;
	//полуциклы
	std::vector<std::shared_ptr<Halfcycle>> halfcycles;
	//интервалы параметра Одквиста
	std::vector<interval> intervals;
	//модуль упругости
	double average_modulus;

	//интерполяционный многочлен лагранжа для первого полуцикла
	double f(double eps);

	//производная полинома
	double fd(double eps);
};


//функция напряжения модели
double sig(double eps, double b, double a, double d, double eps_yield_init, double sig_yield_init, nlopt_copy_utility& s);

class Fatigue {
public:
	Fatigue(std::string filename);
private:
	//полуциклы
	std::vector<std::shared_ptr<Halfcycle>> _halfcycles;
	//интервалы параметра Одквиста
	std::vector<interval> _intervals;
	//кол-во полуциклов
	int _halfcycles_num = 0;
	//усредненный модуль упругости
	double _average_modulus = 0.0;
	
	
	
	//unsigned n, const double *x, double *grad, void* interval_num
	//const std::vector<double> &x, std::vector<double> &grad, void* interval_num
	
};
