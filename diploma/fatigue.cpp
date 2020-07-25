#include "fatigue.h"

Fatigue::Fatigue(std::string filename) {
	int halfcycles_num_whole = 0;
    int halfcycles_num_current = 0;

	_halfcycles.push_back(std::make_shared<Halfcycle>(filename, halfcycles_num_current));
	std::ifstream file(filename);

	do {
	
		double load, strain;
		int num;
		file >> num >> load >> strain;

		if (num > halfcycles_num_current) {
			halfcycles_num_whole++;
			halfcycles_num_current = num;
            bool positive_load = (halfcycles_num_whole % 2 == 0);
			_halfcycles.push_back(std::make_shared<Halfcycle>(filename, halfcycles_num_current,
                _halfcycles.back()->modulus(), positive_load));
		}

	} while (!file.eof());
	file.close();
	halfcycles_num_whole++;
	_halfcycles_num = halfcycles_num_whole;	

	//считаем параметр Одквиста на каждом полуцикле (кроме последнего)
	//_halfcycles[0]->set_odquist(0.0);
	for (int i = 0; i < _halfcycles_num - 1; i++) {
		double eps1 = _halfcycles[i]->zero_stress_deform();
		double eps2 = _halfcycles[i+1]->zero_stress_deform();
		double odqist_gain = abs(eps1 - eps2);
		_halfcycles[i]->set_odquist_gain(odqist_gain);

		if (i == 0)
			_halfcycles[i]->set_odquist(0.0);
		else if (i == 1)
			_halfcycles[i]->set_odquist(odqist_gain + _halfcycles[i - 1]->odquist_gain());
		else
			_halfcycles[i]->set_odquist(odqist_gain + _halfcycles[i - 1]->odquist());

        //std::cout << _halfcycles[i]->odquist() << std::endl;
	}

	//считаем усредненный модуль упругости
	_average_modulus = 0.0;
	for (int i = 0; i < _halfcycles_num; i++)
		_average_modulus += _halfcycles[i]->modulus();
	_average_modulus /= _halfcycles_num;

	//Зададим значения параметра Одквиста через равные промежутки
   /* double h = _halfcycles[_halfcycles.size()-2]->odquist() / 10;

    for (int i = 0; i < 10; i++) {
        _intervals.push_back(interval());
        _intervals[i].odqist_init = h * i;
        _intervals[i].odqist_final = h * (i + 1);

        for (int j = 0; j < _halfcycles_num - 1; j++)
            if ((_halfcycles[j]->odquist() >= _intervals[i].odqist_init) &&
                (_halfcycles[j]->odquist() < _intervals[i].odqist_final))
                _intervals[i].odqist_scale.push_back(j);
    }

	for (int i = 0; i < _halfcycles_num; i++)
		std::cout << "p.s = " << _halfcycles[i]->proportional_stress_point_num() <<
		"; y.s = " << _halfcycles[i]->yield_stress_point_num() << "; of total " <<
		_halfcycles[i]->points_num() << std::endl;*/

    //в логарифмической шкале по основанию 10
    //найдем максимальное значение мараметра Одквиста в шкале
    /*bool flag = false;
    int power = -5;
    int parts_num = 1;
    do {

        bool flag_current = false;

        double current_scale = pow(10.0, power);
        for(int i = 0; i < _halfcycles_num - 1; i++)
            if(_halfcycles[i]->odquist() > current_scale)
                flag_current = true;

        if(!flag_current)
            flag = true;

        if(!flag) {
            power++;
            parts_num++;
        }

    } while (!flag);*/

    /*power = -5;
    _intervals.push_back(interval());
    _intervals[0].odqist_init = 0;
    _intervals[0].odqist_final = pow(10.0, power);

    for (int i = 1; i < parts_num; i++) {
        power++;

		_intervals.push_back(interval());
        _intervals[i].odqist_init = _intervals[i-1].odqist_final;
        _intervals[i].odqist_final = pow(10.0, power);

		for (int j = 0; j < _halfcycles_num - 1; j++)
			if ((_halfcycles[j]->odquist() >= _intervals[i].odqist_init) &&
				(_halfcycles[j]->odquist() < _intervals[i].odqist_final))
                _intervals[i].odqist_scale.push_back(j);
    }*/

    /*_intervals.push_back(interval());
    _intervals[0].odqist_init = 0.01;
    _intervals[0].odqist_final = 0.05;

    _intervals.push_back(interval());
    _intervals[1].odqist_init = 0.05;
    _intervals[1].odqist_final = 0.1;

	_intervals.push_back(interval());
	_intervals[2].odqist_init = 0.1;
	_intervals[2].odqist_final = 0.2;

	_intervals.push_back(interval());
	_intervals[3].odqist_init = 0.2;
	_intervals[3].odqist_final = 0.3;

	_intervals.push_back(interval());
	_intervals[4].odqist_init = 0.3;
	_intervals[4].odqist_final = 0.4;

	_intervals.push_back(interval());
	_intervals[5].odqist_init = 0.4;
	_intervals[5].odqist_final = 0.5;

	_intervals.push_back(interval());
	_intervals[6].odqist_init = 0.5;
	_intervals[6].odqist_final = 1.0;*/

	int inum = _halfcycles_num - 1;
	if(inum % 2 != 0)
		for (int i = 0; i < inum - 1; i += 2) {
			_intervals.push_back(interval());
			_intervals.back().odqist_init = _halfcycles[i]->odquist();
			_intervals.back().odqist_final = _halfcycles[i+2]->odquist();
		}
	else {
		for (int i = 0; i < inum - 2; i += 2) {
			_intervals.push_back(interval());
			_intervals.back().odqist_init = _halfcycles[i]->odquist();
			_intervals.back().odqist_final = _halfcycles[i + 2]->odquist();
		}

		_intervals.push_back(interval());
		_intervals.back().odqist_init = _halfcycles[inum - 2]->odquist();
		_intervals.back().odqist_final = _halfcycles[inum - 1]->odquist();
	}

	std::vector<int> distro;

    for (int i = 0; i < _intervals.size(); i++) {
		distro.push_back(0);
        for (int j = 0; j < _halfcycles_num - 1; j++)
			if ((_halfcycles[j]->odquist() >= _intervals[i].odqist_init) &&
				(_halfcycles[j]->odquist() < _intervals[i].odqist_final)) {
				_intervals[i].odqist_scale.push_back(j);
				distro[i]++;
			}
		std::cout << i << "-th interval: " << distro[i] << " h-cles; odqist: [" << _intervals[i].odqist_init
			<< ", " << _intervals[i].odqist_final << ")" << std::endl;
    }

	std::cout << std::endl;

	//минимизация функционала по каждому интервалу для нахождения параметров a, b на каждом интервале
	//копируем данные
	double MAX_TIME = 10.0;
	nlopt_copy_utility s;
	s.average_modulus = _average_modulus;
	s.intervals = _intervals;
	s.halfcycles = _halfcycles;

	for (int i = 0; i < _intervals.size(); i++) {

		s.current_interval_num = i;
		void *data = &s;
		nlopt::opt opt(nlopt::GN_ESCH, 2);
		//GN_CRS2_LM GN_ESCH
		double half_span = 10.0;
		double points_num = 10000.0;
		std::vector<double> lb = { 0.01, 0.01 };
		std::vector<double> ub = { half_span, half_span };
		std::vector<double> steps = { half_span / points_num, half_span / points_num };
		opt.set_xtol_rel(1e-5);
		opt.set_initial_step(steps);

		opt.set_lower_bounds(lb);
		opt.set_upper_bounds(ub);

		opt.set_min_objective(J, data);
		
		//opt.set_population(1);

		opt.set_maxtime(MAX_TIME);

		double min;
		if (opt.optimize(_intervals[i].params, min) < 0)
			throw std::runtime_error("Неудачная оптимизация");
		else
			std::cout << i << "-th interval: Min = " << min << "; a = " << _intervals[i].params[0] <<
			"; b = " << _intervals[i].params[1] << std::endl;
	}

	std::ofstream myfile_d;
	std::ofstream myfile_a;
	std::ofstream myfile_b;
	std::ofstream myfile_E;
	std::ofstream myfile_yield;
	std::ofstream myfile_odquist;
	myfile_d.open("d.alv");
	myfile_a.open("a.alv");
	myfile_b.open("b.alv");
	myfile_E.open("Eaverage.alv");
	myfile_yield.open("yield.alv");
	myfile_odquist.open("odqist.alv");

	myfile_E << _average_modulus << std::endl;
	
	//myfile_odquist << 0 << std::endl;

	myfile_d << 0.0 << "	" << 1.0 << std::endl;
	for (int i = 0; i < _halfcycles_num - 1; i++) {
		myfile_d << _halfcycles[i]->odquist() << "	" << _halfcycles[i]->d() << std::endl;
		myfile_yield << _halfcycles[i]->yield_deform() << "	" << _halfcycles[i]->yield_stress() << std::endl;
		myfile_odquist << _halfcycles[i]->odquist() << std::endl;
	}

	myfile_a << 0.0 << "	" << 2.0 << std::endl;
	myfile_b << 0.0 << "	" << 2.0 << std::endl;
	for (int i = 0; i < _intervals.size(); i++) {
		myfile_a << (_intervals[i].odqist_final + _intervals[i].odqist_init) / 2 << "	" << _intervals[i].params[0] << std::endl;
		myfile_b << (_intervals[i].odqist_final + _intervals[i].odqist_init) / 2 << "	" << _intervals[i].params[1] << std::endl;
	}
		
	myfile_a.close();
	myfile_b.close();
	myfile_d.close();
	myfile_E.close();
	myfile_yield.close();
	myfile_odquist.close();
    std::cout << "Finished!" << std::endl;
}

double sig(double eps, double a, double b, double d, double eps_yield_init, double sig_yield_init, nlopt_copy_utility& s) {
	double eps_s_star = a / d * eps_yield_init;
	double sign_s_star = a * sig_yield_init;
	double value = 0.0;

	if (eps <= eps_s_star)
		value = s.average_modulus * d * eps;
	else
		value = s.average_modulus * d * eps + d * b * (s.f(eps_yield_init + (eps - eps_s_star) / b) - sig_yield_init);
	return value;
}

double nlopt_copy_utility::f(double eps) {
	int n = halfcycles[0]->points_num();
	double sum = 0.0;
	
	if (eps <= halfcycles[0]->yield_deform())
		sum = halfcycles[0]->modulus() * eps;
	else
		for (int i = 0; i < n; i++) {
			double mult = halfcycles[0]->stress(i);
			for (int j = 0; j < n; j++) {
				if (j == i)
					continue;
				mult *= (eps - halfcycles[0]->deform(j)) / (halfcycles[0]->deform(i) - halfcycles[0]->deform(j));
			}
			sum += mult;
		}
	
	return sum;
}

double nlopt_copy_utility::fd(double eps) {
	int n = halfcycles[0]->points_num();
	double sum = 0.0;

	if (eps <= halfcycles[0]->yield_deform())
		sum = halfcycles[0]->modulus();
	else
		for (int i = 0; i < n; i++) {
			double mult_whole = 0.0;;

			for (int k = 0; k < n; k++) {
				double mult = 1.0;

				for (int j = 0; j < n; j++) {
					if (j == i)
						continue;
					if (j == k)
						mult *= 1.0 / (halfcycles[0]->deform(i) - halfcycles[0]->deform(j));
					else
						mult *= (eps - halfcycles[0]->deform(j)) / (halfcycles[0]->deform(i) - halfcycles[0]->deform(j));
				}

				mult_whole += mult;
			}
			
			mult_whole *= halfcycles[0]->stress(i);
			sum += mult_whole;
		}

	return sum;
}

double J(unsigned n, const double *x, double *grad, void* data) {
	//x[0] - a
	//x[1] - b

	/*if ((abs(x[0]) <= 0.0001) || (abs(x[1]) <= 0.0001))
		return INFINITY;*/

	double a = x[0];
	double b = x[1];

	nlopt_copy_utility s = *(nlopt_copy_utility*)data;
	int i = s.current_interval_num;

	double sum = 0.0;
	double grad_a_sum = 0.0;
	double grad_b_sum = 0.0;

    for (int j = 0; j < s.intervals[i].odqist_scale.size(); j++) {
        int index = s.intervals[i].odqist_scale[j];
		double d = s.halfcycles[index]->d();
		double eps_yield = s.halfcycles[index]->yield_deform();
		double sig_yield = s.halfcycles[index]->yield_stress();
		double eps_yield_init = s.halfcycles[0]->yield_deform();
		double sig_yield_init = s.halfcycles[0]->yield_stress();
		double final_eps = s.halfcycles[index]->final_deform();

		double sig_D = s.halfcycles[index]->final_stress();
		double sig_C = sig(final_eps, x[0], x[1], d, eps_yield_init, sig_yield_init, s);

		double eps_yield_s = x[0] * eps_yield_init / d;
		double sig_yield_s = x[0] * sig_yield_init;

		double eps1 = abs(abs(s.halfcycles[index]->init_deform() - abs(eps_yield)));
		double sig1 = abs(abs(s.halfcycles[index]->init_stress() - abs(sig_yield)));

		//+ или -?
		double part1 = sqrt(abs(eps_yield_s - eps1)*abs(eps_yield_s - eps1) + 
		((sig_yield_s - sig1) / s.average_modulus)*((sig_yield_s - sig1) / s.average_modulus));

		double part2 = abs((sig_C - sig_D) / s.average_modulus);

		sum += (part1 + part2);

		//считаем градиент

		double grad_a1 = (1 / part1) * ((eps_yield_init / (d*d))*(x[0] * eps_yield_init - d * eps1)
			+ (sig_yield_init / (s.average_modulus*s.average_modulus)) * (x[0] * sig_yield_init - sig1));
		double grad_b1 = 0.0;

		double pass_to_poly = eps_yield_init + (final_eps - (x[0] / d)*eps_yield_init) / x[1];

		double grad_a2 = s.average_modulus*eps_yield_init + d * x[1] *
			s.fd(pass_to_poly)*(-eps_yield_init / (d * x[1]));

		double grad_b2 = d * (s.f(pass_to_poly) - sig_yield_init) + d * x[1] * s.fd(pass_to_poly)*
			((x[0] / d)*eps_yield_init - final_eps) / (x[1] * x[1]);

		grad_a_sum += grad_a1 + grad_a2;
		grad_b_sum += grad_b1 + grad_b2;
	}

	if (grad) {
		grad[0] = grad_a_sum;
		grad[1] = grad_b_sum;
	}
	
	/*if(sum < 0.1)
		std::cout << grad_a_sum << " " << grad_b_sum << std::endl;*/

	return sum;
}
