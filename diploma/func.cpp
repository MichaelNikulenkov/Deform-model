#include "func.h"

void deleteEmptyLines(const std::string &FilePath)
{
	std::ifstream in(FilePath);

	if (!in.is_open())
	{
		throw std::runtime_error("Ошибка открытия файла");
	}

	std::string line, text;
	while (std::getline(in, line)) {
		if (!(line.empty() || line.find_first_not_of(' ') == std::string::npos))
			text += line + "\n";
	}
		
	in.close();
	std::ofstream out(FilePath);
	out << text;
	out.close();

	/*std::ifstream fin(FilePath.c_str());

	std::ofstream fout;
	fout.open("Exp/temp.txt", std::ios::out);

	std::string str;
	while (getline(fin, str))
	{
		while (str.length() == 0)
			getline(fin, str);

		fout << str << std::endl;
	}
	fout.close();
	fin.close();
	remove(FilePath.c_str());
	rename("Exp/temp.txt", FilePath.c_str());*/

}

void read_data(dataPoint& original_points, dataPoint& transformed_points, std::string filename) {
	int halfcycles_num_whole = 0;
	int halfcycles_num_current = 0;

	original_points.resize(0);
	transformed_points.resize(0);

	deleteEmptyLines(filename);

	//read all data from file to vectors
	std::vector<int> num_vec;
	std::vector<double> load_vec;
	std::vector<double> strain_vec;

	std::ifstream file(filename);
	do {
		double load, strain;
		int num;
		std::string line;
		std::getline(file, line);
		std::istringstream iss(line);
		iss >> num >> load >> strain;

		//переходим в координаты напряжений и относ. деформаций
		num_vec.push_back(num);
		//load_vec.push_back((load * 1000.0) / 1.4998670186584730623899313606538e-5);
		//strain_vec.push_back((strain / 1000.0) / 0.013 * 100);
		load_vec.push_back(load * 1000);
		strain_vec.push_back((strain /*/ 1.3*/) * 100);
	} while (!file.eof());
	file.close();

	//поправим нумерацию полуциклов
	int current_max_num = 0;
	int k = 0;
	for (int i = 0; i < num_vec.size(); i++) {
		if (num_vec[i] > current_max_num) {
			current_max_num = num_vec[i];
			k++;
			num_vec[i] = k;
		}
		else {
			num_vec[i] = k;
		}
	}

	original_points.resize(k + 1);
	transformed_points.resize(k + 1);

	//заполним данные для полуциклов без преобразования координат
	for (int j = 0; j < original_points.size(); j++)
		for (int i = 0; i < num_vec.size(); i++) {
			if (num_vec[i] == j)
				original_points[j].push_back(Point(strain_vec[i], load_vec[i]));
		}

	//заполним данные для преобразованных координат. Начинаем со второго полуцикла т.к. первый в преобразовании
	//не нуждается
	transformed_points.resize(original_points.size());
	transformed_points[0] = original_points[0];
	for (int i = 1; i < original_points.size(); i++) {
		transformed_points[i].resize(original_points[i].size());
		Point new_cs(original_points[i - 1].back().x, original_points[i - 1].back().y);

		for (int j = 0; j < original_points[i].size(); j++) {
			transformed_points[i][j].x = original_points[i][j].x - new_cs.x;
			transformed_points[i][j].y = original_points[i][j].y - new_cs.y;
			//в каждом втором полуцикле меняется направление осей
			if (i % 2 != 0) {
				transformed_points[i][j].x *= -1.0;
				transformed_points[i][j].y *= -1.0;
			}
		}
	}

	//подкорректируем сдвиги точек
	for (int i = 0; i < transformed_points.size(); i++) {
		if (transformed_points[i][0].x < 0) {
			double shift = abs(transformed_points[i][0].x);
			for (int j = 1; j < transformed_points[i].size(); j++)
				transformed_points[i][j].x += shift;
		}
		else {
			double shift = abs(transformed_points[i][0].x);
			for (int j = 1; j < transformed_points[i].size(); j++)
				transformed_points[i][j].x -= shift;
		}

		if (transformed_points[i][0].y < 0) {
			double shift = abs(transformed_points[i][0].y);
			for (int j = 1; j < transformed_points[i].size(); j++)
				transformed_points[i][j].y += shift;
		}
		else {
			double shift = abs(transformed_points[i][0].y);
			for (int j = 1; j < transformed_points[i].size(); j++)
				transformed_points[i][j].y -= shift;
		}

		transformed_points[i][0].x = 0.0;
		transformed_points[i][0].y = 0.0;
	}


	dataPoint original_points1, transformed_points1;
	original_points1.resize(original_points.size());
	transformed_points1.resize(transformed_points.size());

	for (int i = 0; i < transformed_points.size(); i++)
		for (int j = 1; j < transformed_points[i].size(); j++)
			if (transformed_points[i][j].x > 0 && transformed_points[i][j].y > 0) {
				//transformed_points[i].erase(transformed_points[i].begin() + j);
				//original_points[i].erase(original_points[i].begin() + j);
				original_points1[i].push_back(original_points[i][j]);
				transformed_points1[i].push_back(transformed_points[i][j]);
			}

	original_points.clear();
	transformed_points.clear();

	original_points = original_points1;
	transformed_points = transformed_points1;
}

void read_data_no_scale(dataPoint& original_points, dataPoint& transformed_points, std::string filename) {
	int halfcycles_num_whole = 0;
	int halfcycles_num_current = 0;

	original_points.resize(0);
	transformed_points.resize(0);

	deleteEmptyLines(filename);

	//read all data from file to vectors
	std::vector<int> num_vec;
	std::vector<double> load_vec;
	std::vector<double> strain_vec;

	std::ifstream file(filename);
	do {
		double load, strain;
		int num;
		std::string line;
		std::getline(file, line);
		std::istringstream iss(line);
		iss >> num >> load >> strain;

		//переходим в координаты напряжений и относ. деформаций
		num_vec.push_back(num);
		//load_vec.push_back((load * 1000.0) / 1.4998670186584730623899313606538e-5);
		//strain_vec.push_back((strain / 1000.0) / 0.013 * 100);
		load_vec.push_back(load);
		strain_vec.push_back(strain);
	} while (!file.eof());
	file.close();

	//поправим нумерацию полуциклов
	int current_max_num = 0;
	int k = 0;
	for (int i = 0; i < num_vec.size(); i++) {
		if (num_vec[i] > current_max_num) {
			current_max_num = num_vec[i];
			k++;
			num_vec[i] = k;
		}
		else {
			num_vec[i] = k;
		}
	}

	original_points.resize(k + 1);
	transformed_points.resize(k + 1);

	//заполним данные для полуциклов без преобразования координат
	for (int j = 0; j < original_points.size(); j++)
		for (int i = 0; i < num_vec.size(); i++) {
			if (num_vec[i] == j)
				original_points[j].push_back(Point(strain_vec[i], load_vec[i]));
		}

	//заполним данные для преобразованных координат. Начинаем со второго полуцикла т.к. первый в преобразовании
	//не нуждается
	transformed_points.resize(original_points.size());
	transformed_points[0] = original_points[0];
	for (int i = 1; i < original_points.size(); i++) {
		transformed_points[i].resize(original_points[i].size());
		Point new_cs(original_points[i - 1].back().x, original_points[i - 1].back().y);

		for (int j = 0; j < original_points[i].size(); j++) {
			transformed_points[i][j].x = original_points[i][j].x - new_cs.x;
			transformed_points[i][j].y = original_points[i][j].y - new_cs.y;
			//в каждом втором полуцикле меняется направление осей
			if (i % 2 != 0) {
				transformed_points[i][j].x *= -1.0;
				transformed_points[i][j].y *= -1.0;
			}
		}
	}
}

void calc_current_modulus(const dataPoint& transformed_points, dataDouble& current_modulus) {
	current_modulus.resize(transformed_points.size());
	for (int i = 0; i < transformed_points.size(); i++) {
		current_modulus[i].resize(transformed_points[i].size());
		double modulus_sum_numerator = 0.0;
		double modulus_sum_denominator = 0.0;
		for (int j = 0; j < transformed_points[i].size(); j++) {
			modulus_sum_numerator += transformed_points[i][j].y * transformed_points[i][j].x;
			modulus_sum_denominator += transformed_points[i][j].x * transformed_points[i][j].x;
			current_modulus[i][j] = modulus_sum_numerator / modulus_sum_denominator;
		}
	}
}

void calc_proportional_stress(const dataPoint& transformed_points, const dataDouble& current_modulus, indexData& pr_stress_ind) {

	//вычислим отнотносительная погр-ть между значениями модулей упругости
	dataDouble delta;
	delta.resize(transformed_points.size());
	for (int i = 0; i < transformed_points.size(); i++) {
		delta[i].resize(transformed_points[i].size());
		for (int j = 1; j < transformed_points[i].size(); j++)
			delta[i][j] = abs((current_modulus[i][j] - current_modulus[i][j - 1]) / current_modulus[i][j]);
	}

	//найдем точку начала стабильного участка кривой
	indexData num_A;
	num_A.resize(transformed_points.size());
	for (int i = 0; i < num_A.size(); i++)
		num_A[i] = -1;

	double accuracy = 1e-2;

	
	for (int i = 0; i < transformed_points.size(); i++) {
		for (int j = 3; j < transformed_points[i].size(); j++)
			if (delta[i][j] < accuracy) {
				num_A[i] = j;
				break;
			}
	}

	for(int i = 0; i < num_A.size(); i++)
		if (num_A[i] < 0)
			throw std::runtime_error("Не было найдено начало стабильного участка");
	
	//найдем аппроксимации производных модуля упругости (центральной разностью)
	dataDouble modulus_der;
	modulus_der.resize(transformed_points.size());
	for (int i = 0; i < transformed_points.size(); i++) {
		modulus_der[i].resize(transformed_points[i].size());
		for (int j = 2; j < transformed_points[i].size() - 1; j++)
			modulus_der[i][j] = (current_modulus[i][j+1] - current_modulus[i][j - 1]) / abs(transformed_points[i][j + 1].x - transformed_points[i][j - 1].x);
	}

	//write_to_file(modulus_der[1],"modulus2_der.alv");

	//найдем точку предела пропорциональности (не рассм. последнюю точку тк. считаем центральной разностью)
	//indexData num_B;
	pr_stress_ind.resize(transformed_points.size());
	for (int i = 0; i < transformed_points.size(); i++) {
		pr_stress_ind[i] = num_A[i];

		for (int j = num_A[i] + 1; j < transformed_points[i].size() - 1; j++)
			if (abs(modulus_der[i][pr_stress_ind[i]]) > abs(modulus_der[i][j]))
				pr_stress_ind[i] = j;
	}

	//std::cout << "A " << num_A[0] + 1 << " B " << pr_stress_ind[0] + 1 << std::endl;
}

void calc_modulus(const dataDouble& current_modulus, const indexData& pr_stress_ind, vec& modulus, double& average_modulus) {
	//возьмем модуль упругости из найденного предела пропорциональности
	for (int i = 0; i < current_modulus.size(); i++)
		modulus.push_back(current_modulus[i][pr_stress_ind[i]]);

	double sum = 0.0;
	for (int i = 0; i < modulus.size(); i++)
		sum += modulus[i];
	average_modulus = sum / modulus.size();
}

void calc_odqist(const dataPoint& original_points, vec& odquist) {
	odquist.resize(original_points.size());

	//находим точки пересечения с прямой sigma = 0
	indexData zero_stress_deform_ind;
	zero_stress_deform_ind.resize(original_points.size());

	for (int i = 0; i < original_points.size(); i++) {
		zero_stress_deform_ind[i] = 0;
		for (int j = 1; j < original_points[i].size(); j++)
			if (abs(original_points[i][zero_stress_deform_ind[i]].y) > abs(original_points[i][j].y))
				zero_stress_deform_ind[i] = j;
	}

	//считаем параметр Одквиста на каждом полуцикле (кроме последнего)
	vec odquist_gain;
	odquist_gain.resize(original_points.size());

	for (int i = 0; i < original_points.size(); i++) {
		double eps1, eps2, gain;

		if (!(i == (original_points.size() - 1))) {
			eps1 = original_points[i][zero_stress_deform_ind[i]].x;
			eps2 = original_points[i + 1][zero_stress_deform_ind[i + 1]].x;
			gain = abs(eps1 - eps2);
			odquist_gain[i] = gain;
		}

		if (i == 0)
			odquist[i] = 0.0;
		else
			odquist[i] = odquist[i - 1] + odquist_gain[i - 1];
		
	}

	//на последнем полуцикле
	
}

void calc_yield_stress(const dataPoint& transformed_points, const vec& modulus, const indexData& pr_stress_ind, vec& yield_stress, vec& yield_deform) {
	yield_stress.resize(transformed_points.size());
	yield_deform.resize(transformed_points.size());

	//найдем предел текучести как точку пересечения с прямой y = Ex - shift
	double shift = 0.02;

	for (int i = 0; i < transformed_points.size(); i++) {

		vec h_vec;
		h_vec.resize(transformed_points[i].size());

		//расстояние по оси Y между прямой Ex - shift и точками 
		double h = INFINITY;
		for (int j = pr_stress_ind[i]; j < transformed_points[i].size(); j++) {
			//double a = abs(abs(transformed_points[i][j].y / modulus[i] + shift) - abs(transformed_points[i][j].x));
			//double b = abs(abs(modulus[i] * (transformed_points[i][j].x - shift)) - abs(transformed_points[i][j].y));
			//double h_current = a * b / sqrt(a*a + b * b);
			double y = transformed_points[i][j].y;
			double y_line = (y + shift) / modulus[i];
			h_vec[j] = abs(y - y_line);


			/*if (h_current < h) {
				h = h_current;
				yield_stress_ind[i] = j;
			}*/
		}

		int min_ind = pr_stress_ind[i];
		for (int j = pr_stress_ind[i] + 1; j < transformed_points[i].size(); j++)
			if (h_vec[j] < h_vec[min_ind])
				min_ind = j;

		yield_deform[i] = transformed_points[i][min_ind].x;
		yield_stress[i] = transformed_points[i][min_ind].y;

		/*if ((min_ind + 1) == h_vec.size()) {
			double x1 = transformed_points[i][min_ind - 1].x;
			double y1 = transformed_points[i][min_ind - 1].y;
			double x2 = transformed_points[i][min_ind].x;
			double y2 = transformed_points[i][min_ind].y;

			double a = x1 * (y2 - y1) / (x2 - x1) - y1 - shift;
			double b = (y2 - y1) / (x2 - x1) - modulus[i];
			double x = (a / b);
			yield_deform[i] = x;
			yield_stress[i] = x * modulus[i] - shift;
		} else {
			if ((h_vec[min_ind - 1] > h_vec[min_ind + 1])) {
				double x1 = transformed_points[i][min_ind].x;
				double y1 = transformed_points[i][min_ind].y;
				double x2 = transformed_points[i][min_ind + 1].x;
				double y2 = transformed_points[i][min_ind + 1].y;

				double a = x1 * (y2 - y1) / (x2 - x1) - y1 - shift;
				double b = (y2 - y1) / (x2 - x1) - modulus[i];
				double x = (a / b);
				yield_deform[i] = x;
				yield_stress[i] = x * modulus[i] - shift;
			}
			else {
				double x1 = transformed_points[i][min_ind - 1].x;
				double y1 = transformed_points[i][min_ind - 1].y;
				double x2 = transformed_points[i][min_ind].x;
				double y2 = transformed_points[i][min_ind].y;

				double a = x1 * (y2 - y1) / (x2 - x1) - y1 - shift;
				double b = (y2 - y1) / (x2 - x1) - modulus[i];
				double x = (a / b);
				yield_deform[i] = x;
				yield_stress[i] = x * modulus[i] - shift;

				
			}
		}*/

		if (yield_deform[i] < 0 || yield_stress[i] < 0)
			throw "Ошибка предела текучести";

	}
}

void set_odquist_scale(const vec& odquist, dataDouble& odquist_intervals, dataInt& odquist_scale_ind) {

	int step = 5;
	int last_ind = 0;
	for (int i = 0; (i + step) < odquist.size(); i += step) {
		odquist_intervals.push_back({ 0.0, 0.0 });
		odquist_intervals.back()[0] = odquist[i];
		odquist_intervals.back()[1] = odquist[i + step];

		last_ind = i + step;
	}

	if ((last_ind + 1) != odquist.size()) {
		odquist_intervals.push_back({ 0.0, 0.0 });
		odquist_intervals.back()[0] = odquist[last_ind];
		odquist_intervals.back()[1] = odquist.back();
	}


	//распределим индексы полуциклов по интервалам шкалы 
	odquist_scale_ind.resize(odquist_intervals.size());
	for (int i = 0; i < odquist_intervals.size(); i++)
		for (int j = 1; j < odquist.size(); j++)
			if ((odquist_intervals[i][0] <= odquist[j]) && (odquist[j] < odquist_intervals[i][1]))
				odquist_scale_ind[i].push_back(j);

	odquist_scale_ind.back().push_back(odquist.size() - 1);
}

double f(double deform, const dataPoint& transformed_points, double init_yield_deform, double init_modulus) {
	int n = transformed_points[0].size();
	double sum = 0.0;

	if (deform <= init_yield_deform/*transformed_points[0][yield_stress_ind[0]].x*/)
		sum = init_modulus * deform;
	else if (deform < transformed_points[0].back().x) {
		for (int i = 0; i < n - 1; i++)
			if ((deform >= transformed_points[0][i].x) && (deform < transformed_points[0][i + 1].x)) {
				sum = transformed_points[0][i].y + (deform - transformed_points[0][i].x) * (transformed_points[0][i + 1].y - transformed_points[0][i].y) / (transformed_points[0][i + 1].x - transformed_points[0][i].x);
				break;
			}
	}
	else {
		double x1 = transformed_points[0][transformed_points[0].size() - 2].x;
		double y1 = transformed_points[0][transformed_points[0].size() - 2].y;

		double x2 = transformed_points[0][transformed_points[0].size() - 1].x;
		double y2 = transformed_points[0][transformed_points[0].size() - 1].y;

		sum = (deform - x1) * (y2 - y1) / (x2 - x1) + y1;
		
	}
		

	return sum;
}

double fd(double deform, const dataPoint& transformed_points, double init_yield_deform, double init_modulus) {
	int n = transformed_points[0].size();
	double sum = 0.0;

	if (deform <= init_yield_deform)
		sum = init_modulus;
	else if (deform < transformed_points[0].back().x) {
		for(int i = 0; i < n - 1; i++)
			if ((deform >= transformed_points[0][i].x) && (deform < transformed_points[0][i + 1].x)) {
				sum = (transformed_points[0][i + 1].y - transformed_points[0][i].y)/ (transformed_points[0][i + 1].x - transformed_points[0][i].x);
				break;
			}
	}
	else {
		double x1 = transformed_points[0][transformed_points[0].size() - 2].x;
		double y1 = transformed_points[0][transformed_points[0].size() - 2].y;

		double x2 = transformed_points[0][transformed_points[0].size() - 1].x;
		double y2 = transformed_points[0][transformed_points[0].size() - 1].y;

		sum = (y2 - y1) / (x2 - x1) + y1;
	}
		

	return sum;
}

double model_elastic_stress(double deform, double a, double b, double d, double yield_deform_init,
	double yield_stress_init, double init_modulus, const CopyUtility& s) {

	//double eps_s_star = a / d * yield_deform_init;
	//double sign_s_star = a * yield_stress_init;
	double value = 0.0;

	value = init_modulus * d * deform;

	return value;
}

double model_stress(double deform, double a, double b, double d, double yield_deform_init,
	double yield_stress_init, double init_modulus, const CopyUtility& s) {

	double eps_s_star = a / d * yield_deform_init;
	double sign_s_star = a * yield_stress_init;
	double value = 0.0;

	double pass_to_interp = yield_deform_init + (deform - eps_s_star) / b;

	if (deform <= eps_s_star)
		value = init_modulus * d * deform;
	else
		value = init_modulus * d * eps_s_star + d * b *
		(f(pass_to_interp, s.transformed_points, yield_deform_init, init_modulus) - yield_stress_init);

	return value;
}

double model_stress_alpha(double deform, double a, double b, double d, double yield_deform_init,
	double yield_stress_init, double init_modulus, const CopyUtility& s, double alpha) {

	alpha = 1.0;

	double eps_s_star = a / d * yield_deform_init;
	double sign_s_star = a * yield_stress_init;
	double value = 0.0;

	double pass_to_interp = yield_deform_init + (deform - eps_s_star) / b;

	if (deform <= eps_s_star)
		value = init_modulus * d * deform;
	else
		value = init_modulus * d * eps_s_star + d * b *
		(f(pass_to_interp, s.transformed_points, yield_deform_init, init_modulus) - yield_stress_init);

	return alpha * value;
}

double J(unsigned n, const double *x, double *grad, void* data) {

	//double alpha = 0.05;

	//x[0] - a
	//x[1] - b

	/*if ((abs(x[0]) <= 0.0001) || (abs(x[1]) <= 0.0001))
		return INFINITY;*/

	double a = x[0];
	double b = x[1];

	CopyUtility s = *(CopyUtility*)data;
	int i = s.current_interval_num;

	double sum = 0.0;

	//для проверки в какой части функционала ошибка
	double sum_part1 = 0.0;
	double sum_part2 = 0.0;

	double grad_a_sum = 0.0;
	double grad_b_sum = 0.0;

	//int yield_index_init = s.yield_stress_ind[0];
	for (int j = 0; j < s.odquist_scale_ind[i].size(); j++) {

		int index = s.odquist_scale_ind[i][j];
		double w = 1.0;
		/*if (index % 2 == 0)
			w = 0.2;
		else
			w = 1.0;*/

		//пропускаем полуциклы 
		/*if (index % 2 != 0)
			continue;*/

		//int yield_index = s.yield_stress_ind[index];

		double d = s.d[index];
		double eps_yield_exp = s.yield_deform[index];// s.transformed_points[index][yield_index].x;
		double sig_yield_exp = s.yield_stress[index];//s.transformed_points[index][yield_index].y;
		double eps_yield_init = s.yield_deform[0];//s.transformed_points[0][yield_index_init].x;
		double sig_yield_init = s.yield_stress[0];//s.transformed_points[0][yield_index_init].y;
		double final_eps = s.transformed_points[index].back().x;

		double sig_D = s.transformed_points[index].back().y;
		double sig_C = model_stress(final_eps, x[0], x[1], d, eps_yield_init, sig_yield_init, s.modulus[0], s);

		//модельные деф и напр текучести
		double eps_yield_s = (x[0] / d) * eps_yield_init;
		double sig_yield_s = x[0] * sig_yield_init;
		/*double eps_yield_s = x[0] / d * eps_yield_init;
		double sig_yield_s = s.average_modulus * d * eps_yield_s;*/

		//нужно определить еще и напряжение и деф. в начале ЦИКЛА, для этого условие
		double eps1, sig1;
		/*if ((index % 2 != 0) || index == 0) {
			eps1 = abs(abs(s.original_points[index][0].x - abs(s.original_points[index][yield_index].x)));
			sig1 = abs(abs(s.original_points[index][0].y - abs(s.original_points[index][yield_index].y)));
		}
		else {
			eps1 = abs(abs(s.original_points[index - 1][0].x - abs(s.original_points[index][yield_index].x)));
			sig1 = abs(abs(s.original_points[index - 1][0].y - abs(s.original_points[index][yield_index].y)));
		}	*/
		eps1 = eps_yield_exp;
		sig1 = sig_yield_exp;

		double modulus = s.modulus[0];

		double part1 = sqrt((eps_yield_s - eps1)*(eps_yield_s - eps1) +
			(sig_yield_s - sig1)*(sig_yield_s - sig1) / (modulus*modulus));

		double part2 = abs((w * sig_C - sig_D) / modulus);

		sum_part1 += part1;
		sum_part2 += part2;

		sum += (part1 + part2);

		//считаем градиент

		double grad_a1 = (1 / part1) * ((eps_yield_init / (d*d))*(x[0] * eps_yield_init - d * eps1)
			+ (sig_yield_init / (modulus*modulus)) * (x[0] * sig_yield_init - sig1));
		double grad_b1 = 0.0;


		double pass_to_poly = eps_yield_init + (final_eps - (x[0] / d)*eps_yield_init) / x[1];

		double grad_a2 = w * abs((modulus *eps_yield_init + d * x[1] *
			fd(pass_to_poly, s.transformed_points, s.yield_deform[0], s.modulus[0])*(-eps_yield_init / (d * x[1]))));

		double grad_b2 = w * abs((d * (f(pass_to_poly, s.transformed_points, s.yield_deform[0], s.modulus[0]) - sig_yield_init) + d * x[1] * fd(pass_to_poly, s.transformed_points, s.yield_deform[0], s.modulus[0])*
			((x[0] / d)*eps_yield_init - final_eps) / (x[1] * x[1])));
		
		grad_a_sum += (grad_a1 + grad_a2);
		grad_b_sum += (grad_b1 + grad_b2);
		
	}

	if (grad) {
		grad[0] = grad_a_sum;
		grad[1] = grad_b_sum;
	}

	return sum;
}

void calc_d(const vec& modules, vec& d) {
	d.resize(modules.size());

	d[0] = 1.0;
	for (int i = 1; i < d.size(); i++)
		d[i] = modules[i] / modules[0];
}

void calc_ab(double average_modulus, vec& a, vec& b, CopyUtility& s) {
	//минимизация функционала по каждому интервалу для нахождения параметров a, b на каждом интервале
	//копируем данные
	double MAX_TIME = 5.0;
	
	for (int i = 0; i < s.odquist_intervals.size(); i++) {	

		double alpha_borders_max_inc = 1.0;
		double alpha_borders_growth;

		alpha_borders_growth = (i + 1) * alpha_borders_max_inc / s.odquist_intervals.size();
		

		s.current_interval_num = i;

		void *data = &s;
		nlopt::opt opt(nlopt::GD_STOGO, 2);
		//GN_CRS2_LM GN_ESCH GD_MLSL
		double half_span = 10.0;
		double points_num = 10.0;
		std::vector<double> lb, ub;
		lb = { 0.01, 0.01};
		ub = {half_span, half_span};

		/*if (i % 2 != 0) {
			lb = { 0.01, 0.01, 0.99999 };
			ub = { half_span, half_span, 1.0 + alpha_borders_growth };
		}
			
		else {
			lb = { 0.01, 0.01, 1.0 - alpha_borders_growth };
			ub = { half_span, half_span, 1.000001 };
		}*/
			
		//std::vector<double> steps = { half_span / points_num, half_span / points_num, 0.1 };
		opt.set_xtol_rel(1e-9);
		//opt.set_initial_step(steps);

		opt.set_lower_bounds(lb);
		opt.set_upper_bounds(ub);

		opt.set_min_objective(J, data);

		opt.set_maxtime(MAX_TIME);

		std::vector<double> params = { 1.0, 1.0};
		double min;
		if (opt.optimize(params, min) < 0)
			throw std::runtime_error("Неудачная оптимизация");
		/*else
			std::cout << i << "-th interval: Min = " << min << "; a = " << params[0] <<
			"; b = " << params[1] << "; w = " << params[2] << std::endl;*/

		a.push_back(params[0]);
		b.push_back(params[1]);

		std::cout << i + 1 << " of " << s.odquist_intervals.size() << " intervals" << std::endl;
	}

	double asa;
	asa = 1.0;
}

void clear_file(std::string filename) {
	std::ofstream ofs;
	ofs.open(filename, std::ofstream::out | std::ofstream::trunc);
	ofs.close();
}

void write_to_file(const vec& x, const vec& y, std::string filename) {
	std::ofstream outfile;

	outfile.open(filename/*, std::ios_base::app*/);
	for (int i = 0; i < x.size(); i++)
		outfile << x[i] << "	" << y[i] << std::endl;

	outfile.close();
}

void write_to_file(int ind, const vec& x, const vec& y, std::string filename) {
	std::ofstream outfile;

	outfile.open(filename/*, std::ios_base::app*/);
	for (int i = 0; i < x.size(); i++)
		outfile << ind << "	" << x[i] << "	" << y[i] << std::endl;

	outfile.close();
}

void write_to_file(const vec& x, std::string filename) {
	std::ofstream outfile;

	//outfile << std::setprecision(10);

	outfile.open(filename/*, std::ios_base::app*/);
	for (int i = 0; i < x.size(); i++)
		outfile << x[i]  << std::endl;

	outfile.close();
}

void write_to_file(const dataPoint& x, std::string filename) {
	std::ofstream outfile;

	outfile.open(filename/*, std::ios_base::app*/);
	for(int j = 0; j < x.size(); j++)
		for (int i = 0; i < x[j].size(); i++)
			outfile << j << "	" << x[j][i].x << "	" << x[j][i].y << std::endl;
	

	outfile.close();
}

void write_to_file(const std::vector<Point>& x, std::string filename) {
	std::ofstream outfile;

	outfile.open(filename/*, std::ios_base::app*/);
	for (int j = 0; j < x.size(); j++)
		outfile << x[j].x << "	" << x[j].y << std::endl;


	outfile.close();
}

void write_to_file_reverse(const std::vector<Point>& x, std::string filename) {
	std::ofstream outfile;

	outfile.open(filename/*, std::ios_base::app*/);
	for (int j = 0; j < x.size(); j++)
		outfile << x[j].y << "	" << x[j].x << std::endl;


	outfile.close();
}

void write_to_file_second_coord(const std::vector<Point>& x, std::string filename) {
	std::ofstream outfile;

	outfile.open(filename/*, std::ios_base::app*/);
	for (int j = 0; j < x.size(); j++)
		outfile << x[j].y << std::endl;


	outfile.close();
}

void write_to_file(const dataDouble& x, std::string filename) {
	std::ofstream outfile;

	outfile.open(filename/*, std::ios_base::app*/);
	for(int i = 0; i < x.size(); i++)
		for (int j = 0; j < x[i].size(); j++)
			outfile << (i + 1) << "	" << (j + 1) << "	" << x[i][j] << std::endl;


	outfile.close();
}

void write_to_file(const indexData& x, std::string filename) {
	std::ofstream outfile;

	outfile.open(filename/*, std::ios_base::app*/);
	for (int i = 0; i < x.size(); i++)
		outfile << x[i] << std::endl;

	outfile.close();
}

double interpolate_ab(double x, const vec& param, const vec& odquist_medium) {

	double average_gain = 0;
	for (int i = 0; i < odquist_medium.size(); i++)
		average_gain += odquist_medium[i];
	average_gain /= odquist_medium.size();

	int additional_it = 1;

	//Безье
	int n = odquist_medium.size() - 1 + additional_it;
	int n_old = odquist_medium.size() - 1;
	std::vector<Point> points;
	int steps_num = 500;
	double t_step = 1.0 / steps_num;

	for (int i = 0; i <= steps_num; i++) {
		double t = i * t_step;

		double sum_x = 0.0;
		double sum_y = 0.0;
		for (int k = 0; k <= n; k++) {
			double b = (fact(n) / (fact(k)*fact(n - k))) * pow(t, k) * pow(1 - t, n - k);

			if (k <= n_old) {
				sum_x += odquist_medium[k] * b;
				sum_y += param[k] * b;
			}
			else {
				sum_x += (odquist_medium.back() + 100 * (average_gain) * ((k - n_old) + 1)) * b;
				sum_y += (param.back()) * b;
			}
			
		}
		points.push_back(Point(sum_x, sum_y));
	}
	
	double val = -1.0;
	bool found = false;
	for (int i = 0; i < points.size() - 1; i++)
		if ((x >= points[i].x) && (x < points[i + 1].x)) {
			val = points[i].y + (x - points[i].x) * (points[i + 1].y - points[i].y) / (points[i + 1].x - points[i].x);
			found = true;
			break;
		}

	if (!found) {
		int ind = points.size() - 1;
		//int back = 2;
		//val = points[ind - back].y + (x - points[ind - back].x) * (points[ind].y - points[ind - back].y) / (points[ind].x - points[ind - back].x);
		//val = points[ind].y;
		val = points.back().y;
	}



	return val;
}

double interpolate_d(double x, const vec& d_vec, const vec& odquist, int halcycle_num) {
	//Безье
	
	double average_gain = 0;
	for (int i = 0; i < odquist.size(); i++)
		average_gain += odquist[i];
	average_gain /= odquist.size();

	int additional_it = 1;

	std::vector<Point> points;
	int steps_num = 500;
	double t_step = 1.0 / steps_num;

	vec current_odquist;
	vec current_d;
	if (halcycle_num % 2 == 0) {
		for (int i = 0; i < odquist.size(); i++) 
			if (i % 2 == 0) {
				current_odquist.push_back(odquist[i]);
				current_d.push_back(d_vec[i]);
			}
	}
	else {
		for (int i = 0; i < odquist.size(); i++)
			if (i % 2 != 0) {
				current_odquist.push_back(odquist[i]);
				current_d.push_back(d_vec[i]);
			}
	}

	int n_old = current_odquist.size() - 1;
	int n = current_odquist.size() - 1 + additional_it;

	for (int i = 0; i <= steps_num; i++) {
		double t = i * t_step;
		double sum_x = 0.0;
		double sum_y = 0.0;
		for (int k = 0; k <= n; k++) {
			double b = (fact(n) / (fact(k)*fact(n - k))) * pow(t, k) * pow(1 - t, n - k);
			

			if (k <= n_old) {
				sum_x += current_odquist[k] * b;
				sum_y += current_d[k] * b;
			}
			else {
				sum_x += (current_odquist.back() + 100 * (average_gain)* ((k - n_old) + 1)) * b;
				sum_y += (current_d.back()) * b;
			}
		}

		points.push_back(Point(sum_x, sum_y));
	}

	double val = -1.0;
	bool found = false;
	for (int i = 0; i < points.size() - 1; i++)
		if ((x >= points[i].x) && (x < points[i + 1].x)) {
			val = points[i].y + (x - points[i].x) * (points[i + 1].y - points[i].y) / (points[i + 1].x - points[i].x);
			found = true;
			break;
		}

	if (!found) {
		int ind = points.size()-1;
		//int back = 2;
		//val = points[ind - back].y + (x - points[ind - back].x) * (points[ind].y - points[ind - back].y) / (points[ind].x - points[ind - back].x);
		val = points.back().y;
	}

	return val;
}

long double fact(int N)
{
	if (N < 0) // если пользователь ввел отрицательное число
		return 0; // возвращаем ноль
	if (N == 0) // если пользователь ввел ноль,
		return 1; // возвращаем факториал от нуля - не удивляетесь, но это 1 =)
	else // Во всех остальных случаях
		return N * fact(N - 1); // делаем рекурсию.
}

void calc_model_odqist_given_init_modulus(vec& model_odquist, double E_init, const CopyUtility& s, const vec& a, const vec& b, const vec& alpha, const vec& odquist_medium) {
	
	int additional_iterations = 200;

	model_odquist.resize(s.transformed_points.size() - 1 + additional_iterations);

	//считаем параметр Одквиста на каждом полуцикле (кроме последнего)
	vec odquist_gain;
	odquist_gain.resize(model_odquist.size());

	vec zero_def_elastic_x;
	zero_def_elastic_x.push_back(0.0);

	//int yield_index_init = s.yield_stress_ind[0];
	double eps_yield_init = s.yield_deform[0];//s.transformed_points[0][yield_index_init].x;
	double sig_yield_init = s.yield_stress[0];//s.transformed_points[0][yield_index_init].y;


	//значение деформации и напрядежения в последней точке прошлого полуцикла
	double last_deform = s.transformed_points[0].back().x;
	double last_stress = s.transformed_points[0].back().y;;
	vec last_stress_vec;
	last_stress_vec.resize(model_odquist.size());
	last_stress_vec[0] = last_stress;

	double average_deform = 0.0;
	for (int i = 1; i < s.transformed_points.size(); i++)
		average_deform += s.transformed_points[i].back().x;
	average_deform /= s.transformed_points.size();

	//на первом полуцикле
	double a0 = 2.0;
	double d0 = 1.0;//interpolate_d(s.odquist[0], s.d, s.odquist, 1);
	odquist_gain[0] = -last_stress / (-E_init * d0);
	model_odquist[0] = 0.0;


	
	for (int i = 1; i < s.transformed_points.size() + additional_iterations - 1; i++) {
		double accuracy = 1e-4;

		//одквист
		double eps1, eps2, gain;	

		double sign;
		if (i % 2 != 0)
			sign = -1.0;
		else
			sign = 1.0;


		model_odquist[i] = model_odquist[i - 1] + odquist_gain[i - 1];

		double odquistk = model_odquist[i];
		double ak = interpolate_ab(odquistk, a, odquist_medium);
		double bk = interpolate_ab(odquistk, b, odquist_medium);
		double dk = interpolate_d(odquistk, s.d, s.odquist, i);
		double final_deform;
		/*if (i < s.transformed_points.size())
			final_deform = s.transformed_points[i].back().x;
		else 
			final_deform = average_deform;*/
		final_deform = average_deform;

		//найдем точку пересечения осью абцисс продолжения упругих участков
		//last_stress + sign * E_init * dk * deform = 0;
		zero_def_elastic_x.push_back(-last_stress / (sign * E_init * dk));

		

		//odquist_gain[i] = abs(zero_def_elastic_x[i - 1] - zero_def_elastic_x[i]);


		last_deform += sign * final_deform;
		last_stress += sign * model_stress(final_deform, ak, bk, dk, eps_yield_init, sig_yield_init, E_init, s);
		last_stress_vec[i] = last_stress;

		//найдем приращение деформаций
		double odquist_prev = model_odquist[i - 1];
		double d_prev = interpolate_d(odquist_prev, s.d, s.odquist, i - 1);
		double next_zero_def = -last_stress_vec[i] / (-sign * E_init * d_prev);
		odquist_gain[i] = abs(zero_def_elastic_x.back() - next_zero_def);
	}
}

void calc_model_odqist_given_init_yield_stress(vec& model_odquist, double yield_deform_init, double yield_str_init, const CopyUtility& s, const vec& a, const vec& b, const vec& alpha, const vec& odquist_medium) {

	int additional_iterations = 200;

	double eps_yield_init = yield_deform_init;
	double sig_yield_init = yield_str_init;
	double E_init = s.modulus[0];

	model_odquist.resize(s.transformed_points.size() - 1 + additional_iterations);

	//считаем параметр Одквиста на каждом полуцикле (кроме последнего)
	vec odquist_gain;
	odquist_gain.resize(model_odquist.size());

	vec zero_def_elastic_x;
	zero_def_elastic_x.push_back(0.0);

	//значение деформации и напрядежения в последней точке прошлого полуцикла
	double last_deform = s.transformed_points[0].back().x;
	double last_stress = s.transformed_points[0].back().y;;
	vec last_stress_vec;
	last_stress_vec.resize(model_odquist.size());
	last_stress_vec[0] = last_stress;

	//на первом полуцикле
	double d0 = 1.0;//interpolate_d(s.odquist[0], s.d, s.odquist, 1);
	odquist_gain[0] = -last_stress / (-E_init * d0);
	model_odquist[0] = 0.0;

	double average_deform = 0.0;
	for (int i = 1; i < s.transformed_points.size(); i++)
		average_deform += s.transformed_points[i].back().x;
	average_deform /= s.transformed_points.size();


	for (int i = 1; i < s.transformed_points.size() + additional_iterations - 1; i++) {
		double accuracy = 1e-4;

		//одквист
		double eps1, eps2, gain;

		double sign;
		if (i % 2 != 0)
			sign = -1.0;
		else
			sign = 1.0;

		model_odquist[i] = model_odquist[i - 1] + odquist_gain[i - 1];

		double odquistk = model_odquist[i];
		double ak = interpolate_ab(odquistk, a, odquist_medium);
		double bk = interpolate_ab(odquistk, b, odquist_medium);
		double dk = interpolate_d(odquistk, s.d, s.odquist, i);
		double final_deform;
		/*if (i < s.transformed_points.size())
			final_deform = s.transformed_points[i].back().x;
		else
			final_deform = average_deform;*/
		final_deform = average_deform;

		//найдем точку пересечения осью абцисс продолжения упругих участков
		//last_stress + sign * E_init * dk * deform = 0;
		zero_def_elastic_x.push_back(-last_stress / (sign * E_init * dk));

		//odquist_gain[i] = abs(zero_def_elastic_x[i - 1] - zero_def_elastic_x[i]);


		last_deform += sign * final_deform;
		last_stress += sign * model_stress(final_deform, ak, bk, dk, eps_yield_init, sig_yield_init, E_init, s);
		last_stress_vec[i] = last_stress;

		//найдем приращение деформаций
		double odquist_prev = model_odquist[i - 1];
		double d_prev = interpolate_d(odquist_prev, s.d, s.odquist, i - 1);
		double next_zero_def = -last_stress_vec[i] / (-sign * E_init * d_prev);
		odquist_gain[i] = abs(zero_def_elastic_x.back() - next_zero_def);

	}

	//проверка
	/*dataPoint pnts;
	pnts.resize(model_odquist.size());
	int size = 100;
	
	double last_deform1 = 0.0;
	double last_stress1 = 0.0;
	for (int i = 0; i < model_odquist.size(); i++) {

		double sign;
		if (i % 2 != 0)
			sign = -1.0;
		else
			sign = 1.0;

		double final_deform;
		if (i < s.transformed_points.size())
			final_deform = s.transformed_points[i].back().x;
		else
			final_deform = average_deform;

		double step = final_deform / size;

		for (int j = 0; j < size; j++) {
			double odquistk = model_odquist[i];
			double ak = interpolate_ab(odquistk, a, odquist_medium);
			double bk = interpolate_ab(odquistk, b, odquist_medium);
			double alphak = interpolate_ab(odquistk, alpha, odquist_medium);
			double dk = interpolate_d(odquistk, s.d, s.odquist, i);

			double def = sign * j * step + last_deform1;
			double val = sign * model_stress_alpha(j * step, ak, bk, dk, eps_yield_init, sig_yield_init, E_init, s, alphak) + last_stress1;
			pnts[i].push_back(Point(def, val));
		}

		double odquistk = model_odquist.back();
		double ak = interpolate_ab(odquistk, a, odquist_medium);
		double bk = interpolate_ab(odquistk, b, odquist_medium);
		double alphak = interpolate_ab(odquistk, alpha, odquist_medium);
		double dk = interpolate_d(odquistk, s.d, s.odquist, i);
		last_deform1 += sign * final_deform;
		last_stress1 += sign * model_stress_alpha(final_deform, ak, bk, dk, eps_yield_init, sig_yield_init, E_init, s, alphak);
	}

	write_to_file(pnts, "modelled_with_odquist.alv");*/
}

double calc_model_prop_and_yield_stress(double E_init, const CopyUtility& s) {

	//int yield_index_init = s.yield_stress_ind[0];
	//double eps_yield_init = s.transformed_points[0][yield_index_init].x;
	//double sig_yield_init = s.transformed_points[0][yield_index_init].y;

	////промоделируем кривые с шагом
	//
	//dataPoint points;
	//points.resize(s.transformed_points.size());
	//
	//for (int i = 0; i < points.size(); i++) {
	//	
	//	double odquistk = s.odquist[i];
	//	double ak = interpolate_ab(odquistk, a, odquist_medium);
	//	double bk = interpolate_ab(odquistk, b, odquist_medium);
	//	double dk = interpolate_d(odquistk, s.d, s.odquist, i);

	//	for (int j = 0; j < s.transformed_points[i].size(); j++) {
	//		double deform = s.transformed_points[i][j].x;
	//		points[i].push_back(Point(deform, model_stress(deform, ak, bk, dk, eps_yield_init, sig_yield_init, E_init, s)));
	//	}
	//}


	dataDouble current_modulus;
	indexData pr_stress_ind;
	vec yield_stress, yield_deform;
	vec modulus;
	double average_modulus;

	calc_current_modulus(s.transformed_points, current_modulus);
	calc_proportional_stress(s.transformed_points, current_modulus, pr_stress_ind);
	calc_modulus(current_modulus, pr_stress_ind, modulus, average_modulus);
	modulus[0] = E_init;
	calc_yield_stress(s.transformed_points, modulus, pr_stress_ind, yield_stress, yield_deform);

	return yield_stress[0];//s.transformed_points[0][yield_stress_ind[0]].y;

	/*pr_stress.resize(points.size());
	yield_stress.resize(points.size());

	for (int i = 0; i < points.size(); i++) {
		pr_stress[i] = points[i][pr_stress_ind[i]].y;
		yield_stress[i] = points[i][yield_stress_ind[i]].y;
	}*/
}


double interpolate_w(const vec& w, const dataInt& odquist_scale_ind, int halfcycle_ind) {

	bool found = false;
	double res;

	//в полуциклах сжатия веса = 1
	if ((halfcycle_ind % 2 != 0) || (halfcycle_ind == 0)) {
		found = true;
		res = 1.0;
	}

	if(!found)
		for (int i = 0; i < odquist_scale_ind.size(); i++) {
			for (int j = 0; j < odquist_scale_ind[i].size(); j++) {
				if (halfcycle_ind == odquist_scale_ind[i][j]) {
					res = w[i];
					found = true;
				}

				if (found)
					break;
			}

			if (found)
				break;
		}
	
	return res;
}

double destruct_crit(double odqist, double gamma, double delta) {
	return pow(odqist / delta, gamma);
}

void calc_with_destruction_crit(indexData& destruct_ind_distro, dataDouble& model_odquist_destruct, const dataDouble& model_odquist) {
	/*double gamma = 1.19614;
	double delta = 0.728882;*/
	double gamma = 1.452;
	double delta = 0.785;

	for (int i = 0; i < model_odquist.size(); i++) {

		bool finished = false;
		int intersect_ind;
		for (int j = 1; j < model_odquist[i].size() - 1; j++) {

			double x1 = model_odquist[i][j];
			double x2 = model_odquist[i][j + 1];
			double y1 = j + 1;
			double y2 = j + 2;
			double y1m = destruct_crit(x1, gamma, delta) + 1.0;
			double y2m = destruct_crit(x2, gamma, delta) + 1.0;

			//double intersect_y = (y1*(y2m - y1m) - y1m * (y2 - y1)) / ((y2m - y2) - (y1m - y1));

			//int y1 = (int)std::floor(y);
			//int y2 = (int)std::ceil(y);
			//std::cout << "x = " << x << " " << "y = " << y << " " << y1 << " " << y2 << std::endl;
			//std::cout << intersect_y << std::endl;

			bool first = (((y1m - y1) > 0) && ((y2m - y2) < 0));
			bool second = (((y1m - y1) < 0) && ((y2m - y2) > 0));

			if (first || second) {
				finished = true;
				if(abs(y1m - y1) < abs(y2m - y2))
					intersect_ind = j;
				else
					intersect_ind = j + 1;
			} 

			if (finished)
				break;
		}

		if (finished) {
			vec v;
			for (int k = 0; k <= intersect_ind; k++)
				v.push_back(model_odquist[i][k]);
			model_odquist_destruct.push_back(v);
			destruct_ind_distro.push_back(intersect_ind);
		}
		else {
			model_odquist_destruct.push_back(model_odquist[i]);
		}

	}
}

void calc_with_destruction_crit_continuous(vec& destruct_ind_distro, dataDouble& model_odquist_destruct, const dataDouble& model_odquist) {
	double gamma = 1.452;
	double delta = 0.785;

	for (int i = 0; i < model_odquist.size(); i++) {

		bool finished = false;
		double intersect_ind;
		for (int j = 1; j < model_odquist[i].size() - 1; j++) {

			double y1 = model_odquist[i][j];
			double y2 = model_odquist[i][j + 1];
			double x1 = j + 1;
			double x2 = j + 2;
			double x1m = destruct_crit(y1, gamma, delta) + 1.0;
			double x2m = destruct_crit(y2, gamma, delta) + 1.0;
			double y1m = y1;
			double y2m = y2;

			//double intersect_y = (y1*(y2m - y1m) - y1m * (y2 - y1)) / ((y2m - y2) - (y1m - y1));

			//int y1 = (int)std::floor(y);
			//int y2 = (int)std::ceil(y);
			//std::cout << "x = " << x << " " << "y = " << y << " " << y1 << " " << y2 << std::endl;
			//std::cout << intersect_y << std::endl;

			bool first = (((x1m - x1) >= 0) && ((x2m - x2) <= 0));
			bool second = (((x1m - x1) < 0) && ((x2m - x2) > 0));

			bool flag = ((x1m <= x1) && (x2m >= x2)) || ((x1m >= x1) && (x2m <= x2));

			if (flag/*first || second*/) {
				finished = true;
				/*if (abs(y1m - y1) < abs(y2m - y2))
					intersect_ind = j;
				else
					intersect_ind = j + 1;*/

				//метод вилки
				double accuracy = 1e-4;
				bool stop = false;

				double x_left = x1m;
				double x_right = x2m;
				
				while (!stop) {
					double x_mid = (x_left + x_right) * 0.5;
					double d = (y2 - y1) / (x2 - x1);
					double val_left = delta * pow(x_left, 1.0 / gamma) - (y1 + d * (x_left - x1));
					double val_mid = delta * pow(x_mid, 1.0 / gamma) - (y1 + d * (x_mid - x1));

					if (val_left*val_mid < 0)
						x_right = x_mid;
					else
						x_left = x_mid;

					if (abs(x_left - x_right) < accuracy)
						stop = true;
				}

				intersect_ind = (x_left + x_right)*0.5;
				//std::cout << intersect_ind << std::endl;
			}

			if (finished)
				break;
		}

		if (finished) {
			vec v;
			for (int k = 0; k <= floor(intersect_ind); k++)
				v.push_back(model_odquist[i][k]);

			/*double intersect_odq = destruct_crit(intersect_ind, gamma, delta) + 1.0;
			

			vec last_point;
			last_point.push_back(intersect_ind);
			last_point.push_back(intersect_odq);
			v.push_back(Point());*/

			model_odquist_destruct.push_back(v);
			destruct_ind_distro.push_back(intersect_ind);
		}
		else {
			model_odquist_destruct.push_back(model_odquist[i]);
		}

	}
}

double deformations_range(dataPoint& original_points) {

	indexData zero_stress_deform_ind;
	zero_stress_deform_ind.resize(original_points.size());

	for (int i = 0; i < original_points.size(); i++) {
		zero_stress_deform_ind[i] = 0;
		for (int j = 1; j < original_points[i].size(); j++)
			if (abs(original_points[i][zero_stress_deform_ind[i]].y) > abs(original_points[i][j].y))
				zero_stress_deform_ind[i] = j;
	}

	vec zero_stress;
	for (int i = 0; i < original_points.size(); i++)
		zero_stress.push_back(original_points[i][zero_stress_deform_ind[i]].x);

	double max_deform = -1.0;
	for (int i = 0; i < original_points.size(); i += 2) {
		double a = abs(zero_stress[i] - zero_stress[i + 1]);
		if (max_deform < a)
			max_deform = a;
	}
		
	return max_deform;
}

double interpolate(double x, vecPoint points) {
	double val = -1.0;
	bool found = false;
	for (int i = 0; i < points.size() - 1; i++)
		if ((x >= points[i].x) && (x < points[i + 1].x)) {
			val = points[i].y + (x - points[i].x) * (points[i + 1].y - points[i].y) / (points[i + 1].x - points[i].x);
			found = true;
			break;
		}

	if (!found) {
		if(x >= points.back().y)
			val = points.back().y;
		else
			val = points[0].y;
	}

	return val;
}

Point calc_criteria_point(std::string name1) {
	dataPoint original_points;
	dataPoint transformed_points;
	vecPoint criteria_point;
	indexData yield_stress_ind;
	dataInt odquist_scale_ind;
	dataDouble odquist_intervals;
	vec odquist;
	vec d;
	vec modulus;
	dataDouble current_modulus;
	indexData pr_stress_ind;
	double average_modulus;
	int current_interval_num;

	std::string name = std::string("C:/Users/we/source/repos/diploma/diploma/Exp/") + name1 +std::string(".alv");
	read_data(original_points, transformed_points, name);
	calc_current_modulus(transformed_points, current_modulus);
	calc_proportional_stress(transformed_points, current_modulus, pr_stress_ind);
	calc_modulus(current_modulus, pr_stress_ind, modulus, average_modulus);
	calc_d(modulus, d);

	calc_odqist(original_points, odquist);

	return Point(log(transformed_points.size()), log(odquist.back()));
}