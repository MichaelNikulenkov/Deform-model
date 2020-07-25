#pragma once
#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include "math.h"

struct point {
	double eps = 0.0; //отн. деформация
	double sig = 0.0; //напряжение
	double curr_modulus = 0.0; //модуль упругости в точке.
	double delta; //отн. погрешность между значениями модуля упругости в тек. и прошлой точке
	double der; //аппроксимация значения производной модуля упругости в точке
};

class Halfcycle {
public:
	//считать полуцикл под номером num и посчитать для него модуль упругости. 
	//init_modulus - мод. упругости на первом полуцикле.
	//x_orig, y _orig - основание локальной системы координат полуцикла
	//inverted_axis - определяет направления осей локальной системы координат относительно глобальной с.к.
    Halfcycle(std::string filename, unsigned int num, double x_orig, double y_orig, bool inverted_axis = false, double init_modulus = -1.0, bool positive_load = true);

	//get-методы
	double deform(int index);
	double stress(int index);
	int proportional_stress_point_num();
	int yield_stress_point_num();
	int points_num();
	double odquist();
	double odquist_gain();
	double modulus();
	double zero_stress_deform();
	double d();
	double yield_stress();
	double yield_deform();
	double proportional_stress();
	double init_stress();
	double init_deform();
	double final_stress();
	double final_deform();

	//set-методы
	void set_odquist(double value);
	void set_odquist_gain(double value);
private:
	//точки
	std::vector<point> _points;
	std::vector<point> _transformed_points;
	//номер точки с пределом пропорциональности
	int _proportional_stress_point_num = 0;
	//точка с нулевым напряжением
	int _zero_stress_deform_point_num = 0;
	//номер точки с пределом текучести
	int _yield_stress_point_num = 0;
	//параметр Одквиста
	double _odquist = 0.0;
	//прирост параметра Одквиста
	double _odquist_gain = 0.0;
	//деформация при нулевом напряжении
	double _zero_stress_deform = 0.0;
	//предел текучести
	double _yield_stress = 0.0;
	//деформация в пределе текучести
	double _yield_deform = 0.0;
	//модуль упругости
	double _modulus = 0.0;
	//остальные параметры полуцикла
	double _d = 0.0;
	//локальная система координат полуцикла
	double _x_origin = 0.0;
	double _y_origin = 0.0;
	//расчет параметров
	void _calc_proportional_stress_point_num();
	void _calc_zero_stress_deform();
	void _calc_a(double _average_modulus);
};
