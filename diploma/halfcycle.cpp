#include "halfcycle.h"

Halfcycle::Halfcycle(std::string filename, unsigned int num, double x_orig, double y_orig, bool inverted_axis = false, double init_modulus = -1.0, bool positive_load = true) {
	if (num > 0 && (init_modulus < 0))
		throw std::runtime_error("Некорректный параметр: начальный модуль упругости");

	double modulus_sum_numerator = 0.0;
	double modulus_sum_denominator = 0.0;
    unsigned int halfcylcle_num = 0;
	std::ifstream file(filename);

	file.clear();
	file.seekg(0, std::ios::beg);

	//считываем точки из файла и считаем параметры для них
	do {
        double load = 0;
        double strain = 0;
		file >> halfcylcle_num >> load >> strain;

		if (halfcylcle_num == num) {

			point v;
			//переходим в координаты напряжений и относ. деформаций
            v.sig = (load * 1000) / 1.4998670186584730623899313606538e-5;
			v.eps = (strain / 1000) / 0.013 * 100;

			//преобразуем к локальным координатам
			point v_loc;
			double a_x = v.eps;
			double a_y = v.sig;

			double c_x = x_orig;
			double c_y = y_orig;

			v_loc.sig = a_x - c_x;
			v_loc.eps = a_y - c_y;

			if (inverted_axis) {
				v_loc.sig *= -1.0;
				v_loc.eps *= -1.0;
			}

			//расчет модуля упругой деформации для каждой точки в виде отношения сумм
			modulus_sum_numerator += v.sig * v.eps;
			modulus_sum_denominator += v.eps*v.eps;
			v.curr_modulus = modulus_sum_numerator / modulus_sum_denominator;
			
			if (_points.size() > 0) {
				//отнотносительная погр-ть между значениями модулей упругости на тек. и предыдущ шагах
				v.delta = abs((v.curr_modulus - _points.back().curr_modulus) / v.curr_modulus);	
			}
			
			_points.push_back(v);
		}	
	} while (halfcylcle_num <= num && !file.eof());

	_calc_zero_stress_deform();
	_calc_proportional_stress_point_num();
	//возьмем модуль упругости из найденного предела пропорциональности
	_modulus = _points[_proportional_stress_point_num].curr_modulus;

	//найдем d
	if (num == 0)
		_d = 1.0;
	else
		_d = _modulus / init_modulus;

    //найдем предел текучести как ближайшую точку к прямой y = Ex -+ shift
	double h = INFINITY;
	double shift = abs(0.002 - abs(_zero_stress_deform));
	int nearest_num = 0;
    double a = 0;
    double b = 0;

	bool negative_start = (_zero_stress_deform < 0);
	if (num == 0)
		negative_start = false;

	for (int i = _proportional_stress_point_num; i < _points.size(); i++) {
        if (!negative_start) {
			a = abs(abs(_points[i].sig/_modulus + shift) - abs(_points[i].eps));
			b = abs(abs(_modulus * (_points[i].eps - shift)) - abs(_points[i].sig));
        } else {
			a = abs(abs(_points[i].sig / _modulus - shift) - abs(_points[i].eps));
			b = abs(abs(_modulus * (_points[i].eps + shift)) - abs(_points[i].sig));
        }

		double h_current = a * b / sqrt(a*a + b*b);
		if (h_current < h) {
			h = h_current;
			nearest_num = i;
			//std::cout << h_current << std::endl;
		}
	}

	_yield_stress_point_num = nearest_num;
	_yield_stress = _points[_yield_stress_point_num].sig;
	_yield_deform = _points[_yield_stress_point_num].eps;

	std::ofstream myfile;
	myfile.open("E.alv");
	for (int i = 0; i < _points.size(); i++)
		myfile << _points[i].eps << "	" << _points[i].curr_modulus << std::endl;
	myfile.close();

	/*std::cout << "E: "<< _modulus << std::endl;
	std::cout <<"Предел пропорциональности: " << "num: " << _proportional_stress_point_num << "; " << _points[_proportional_stress_point_num].eps <<
		", " << _points[_proportional_stress_point_num].sig << std::endl;*/

	file.close();

}

void Halfcycle::_calc_zero_stress_deform() {
	//находим накопленную деформацию как точку пересечения с прямой sigma = 0
	int x_intersect_point = 0;
	for (int i = 1; i < _points.size(); i++)
		if (abs(_points[x_intersect_point].sig) > abs(_points[i].sig))
			x_intersect_point = i;

	_zero_stress_deform = _points[x_intersect_point].eps;
	_zero_stress_deform_point_num = x_intersect_point + 1;
}

void Halfcycle::_calc_proportional_stress_point_num() {
	//найдем точку начала стабильного участка кривой
	int num_A = -1;
    double accuracy = 1e-3;

	for (int i = _zero_stress_deform_point_num - 1; i < _points.size(); i++)
		if (_points[i].delta < accuracy) {
			num_A = i;
			break;
		}

	if (num_A < 0)
		throw std::runtime_error("Не было найдено начало стабильного участка");

	//найдем аппроксимации производных модуля упругости (центральной разностью)
	for (int i = 1; i < _points.size() - 1; i++) {
		_points[i].der = (_points[i + 1].curr_modulus - _points[i - 1].curr_modulus) / abs(_points[i + 1].eps - _points[i - 1].eps);
		/*std::cout << "Е: " << _points[i].curr_modulus << " "
			<< "delta:" << _points[i].delta << " " << "der:" << _points[i].der << std::endl;*/
	}

	//найдем точку предела пропорциональности (не рассм. последнюю точку тк. считаем центральной разностью)
	int num_B = num_A;
	for (int i = num_A + 1; i < _points.size() - 1; i++)
		if (abs(_points[num_B].der) > abs(_points[i].der))
			num_B = i;

	/*std::cout << "Конец осциляций: " << "num: " << num_A << "; " << _points[num_A].eps <<
		", " << _points[num_A].sig << std::endl;*/

	_proportional_stress_point_num = num_B;
}

void Halfcycle::_calc_a(double _average_modulus) {

}

double Halfcycle::modulus() {
	return _modulus;
}

double Halfcycle::d() {
	return _d;
}

double Halfcycle::odquist() {
	return _odquist;
}

double Halfcycle::odquist_gain() {
	return _odquist_gain;
}

void Halfcycle::set_odquist(double val) {
	_odquist = val;
}

void Halfcycle::set_odquist_gain(double val) {
	_odquist_gain = val;
}

double Halfcycle::zero_stress_deform() {
	return _zero_stress_deform;
}

double Halfcycle::init_stress() {
	return _points[0].sig;
}

double Halfcycle::init_deform() {
	return _points[0].eps;
}


double Halfcycle::deform(int index) {
	return _points[index].eps;
}

double Halfcycle::stress(int index) {
	return _points[index].sig;
}

int Halfcycle::points_num() {
	return _points.size();
}

double Halfcycle::yield_stress() {
	return _yield_stress;
}

double Halfcycle::yield_deform() {
	return _yield_deform;
}

double Halfcycle::final_stress() {
	return _points.back().sig;
}

double Halfcycle::final_deform() {
	return _points.back().eps;
}

double Halfcycle::proportional_stress() {
	return _points[_proportional_stress_point_num].sig;
}

int Halfcycle::proportional_stress_point_num() {
	return _proportional_stress_point_num;
}

int Halfcycle::yield_stress_point_num() {
	return _yield_stress_point_num;
}