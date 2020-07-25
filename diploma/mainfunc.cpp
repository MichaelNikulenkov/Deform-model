#include "mainfunc.h"

void calc_everything(const char* FilePath, CopyUtility& s, vec& a, vec& b, vec& alpha, vec& odquist_medium) {
	try {
		dataPoint original_points;
		dataPoint transformed_points;
		vecPoint criteria_point;
		//indexData yield_stress_ind;
		vec yield_stress, yield_deform;
		dataInt odquist_scale_ind;
		dataDouble odquist_intervals;
		vec odquist;
		vec d;
		vec modulus;
		dataDouble current_modulus;
		indexData pr_stress_ind;
		double average_modulus;
		int current_interval_num;

		std::string root = FilePath + std::string("/points");
		_mkdir(root.c_str());

		std::string name = FilePath + std::string(".alv");
		read_data(original_points, transformed_points, name);
		calc_current_modulus(transformed_points, current_modulus);
		calc_proportional_stress(transformed_points, current_modulus, pr_stress_ind);
		calc_modulus(current_modulus, pr_stress_ind, modulus, average_modulus);
		calc_d(modulus, d);

		calc_odqist(original_points, odquist);
		calc_yield_stress(transformed_points, modulus, pr_stress_ind, yield_stress, yield_deform);
		set_odquist_scale(odquist, odquist_intervals, odquist_scale_ind);

		criteria_point.push_back(Point(log(transformed_points.size()-1), log(odquist.back())));
		
		write_to_file(criteria_point, FilePath + std::string("/points/criteria_point.alv"));

		std::cout << name << std::endl;
		std::cout << "Number of halfcycles: " << transformed_points.size() << std::endl;
		std::cout << "Init data: E = " << modulus[0] << "; yield stress = " << yield_stress[0] << std::endl;

		s.transformed_points = transformed_points;
		s.original_points = original_points;
		s.yield_stress = yield_stress;
		s.yield_deform = yield_deform;
		s.odquist_scale_ind = odquist_scale_ind;
		s.odquist_intervals = odquist_intervals;
		s.odquist = odquist;
		s.d = d;
		s.modulus = modulus;
		s.average_modulus = average_modulus;
		s.current_interval_num = 0;

		calc_ab(average_modulus, a, b, s);

		double eps_yield_init = s.yield_deform[0];
		double sig_yield_init = s.yield_stress[0];

		//посчитаем усредненый параметр одквиста по интервалам в шкале (возьмем значение на правой границе)
		odquist_medium.resize(odquist_intervals.size());
		for (int i = 0; i < odquist_intervals.size(); i++)
			odquist_medium[i] = odquist_intervals[i][1];


		//добавим начальные значени€ параметров a и b 
		a.insert(a.begin(), 2.0);
		b.insert(b.begin(), 2.0);
		alpha.insert(alpha.begin(), 1.0);
		odquist_medium.insert(odquist_medium.begin(), 0.0);

		//запишем в файлы параметры
		write_to_file(odquist_medium, a, FilePath + std::string("/points/a.alv"));
		write_to_file(odquist_medium, b, FilePath + std::string("/points/b.alv"));
		write_to_file(odquist, d, FilePath + std::string("/points/d.alv"));
		write_to_file(odquist_medium, FilePath + std::string("/points/odquist_medium.alv"));
		write_to_file(transformed_points, FilePath + std::string("/points/transformed_points.alv"));

		write_to_file(current_modulus[1], FilePath + std::string("/points/E1.alv"));

		//моделирование
		std::vector<vec> diff; //разность между моделью и экспериментом в точках
		dataPoint model_points;
		diff.resize(transformed_points.size());
		model_points.resize(transformed_points.size());
		for (int j = 0; j < transformed_points.size(); j++) {
			int k = j; //номер моделируемого полуцикла 0,1, ...
			double odquistk = odquist[k];
			double ak = interpolate_ab(odquistk, a, odquist_medium);
			double bk = interpolate_ab(odquistk, b, odquist_medium);
			double dk = interpolate_d(odquistk, d, odquist, k);

			//посчитаем разность между моделью и экспериментом в точках
			for (int i = 0; i < transformed_points[k].size(); i++) {
				
				double deform = transformed_points[k][i].x; //ƒеформаци€
				double model = model_stress(deform, ak, bk, dk, eps_yield_init, sig_yield_init, modulus[0], s);

				diff[j].push_back(abs(model - transformed_points[k][i].y));
				model_points[j].push_back(Point(deform, model));
			}
		}


		//запишем в файл
		write_to_file(model_points, FilePath + std::string("/points/model_points.alv"));

		//проверка интерпол€ций
		std::vector<Point> interp_f, interp_fd, interp_a, interp_b, interp_d1, interp_d2;

		vec odquist1, odquist2;
		vec d_vec1, d_vec2;
		for (int i = 0; i < odquist.size(); i++)
			if (i % 2 == 0) {
				odquist1.push_back(odquist[i]);
				d_vec1.push_back(d[i]);
			}
		for (int i = 0; i < odquist.size(); i++)
			if (i % 2 != 0) {
				odquist2.push_back(odquist[i]);
				d_vec2.push_back(d[i]);
			}

		double h_f = transformed_points[0].back().x / (transformed_points[0].size() * 2.0);
		for (int i = 0; i < (transformed_points[0].size() * 2 + 10); i++) {
			interp_f.push_back(Point(h_f*i, f(h_f*i, s.transformed_points, s.yield_deform[0], s.modulus[0])));
			interp_fd.push_back(Point(h_f*i, fd(h_f*i, s.transformed_points, s.yield_deform[0], s.modulus[0])));
		}

		double h_ab = odquist_medium.back() / (odquist.size() * 2.0);
		for (int i = 0; i < ((odquist.size() * 2) + 100); i++) {
			interp_a.push_back(Point(h_ab*i, interpolate_ab(h_ab*i, a, odquist_medium)));
			interp_b.push_back(Point(h_ab*i, interpolate_ab(h_ab*i, b, odquist_medium)));
			interp_d1.push_back(Point(h_ab*i, interpolate_ab(h_ab*i, d_vec1, odquist1)));
			interp_d2.push_back(Point(h_ab*i, interpolate_ab(h_ab*i, d_vec2, odquist2)));
		}

		write_to_file(interp_f, FilePath + std::string("/points/interp_f.alv"));
		write_to_file(interp_fd, FilePath + std::string("/points/interp_fd.alv"));
		write_to_file(interp_a, FilePath + std::string("/points/interp_a.alv"));
		write_to_file(interp_b, FilePath + std::string("/points/interp_b.alv"));
		write_to_file(interp_d1, FilePath + std::string("/points/interp_d1.alv"));
		write_to_file(interp_d2, FilePath + std::string("/points/interp_d2.alv"));
		write_to_file(odquist, FilePath + std::string("/points/odquist.alv"));
		write_to_file(original_points, FilePath + std::string("/points/original_points.alv"));


		//тест преобразовани€
		dataPoint points;
		points.resize(s.transformed_points.size());
		double last_deform1 = s.transformed_points[0].back().x;
		double last_stress1 = s.transformed_points[0].back().y;;
		points[0] = s.transformed_points[0];

		for (int i = 1; i < s.transformed_points.size(); i++) {

			double sign;
			if (i % 2 != 0)
				sign = -1.0;
			else
				sign = 1.0;

			double odquistk = s.odquist[i];
			double ak = interpolate_ab(odquistk, a, odquist_medium);
			double bk = interpolate_ab(odquistk, b, odquist_medium);
			double dk = interpolate_d(odquistk, s.d, s.odquist, i);

			for (int j = 0; j < s.transformed_points[i].size(); j++) {
				double xval = sign * s.transformed_points[i][j].x + last_deform1;
				double yval = sign * model_stress(s.transformed_points[i][j].x, ak, bk, dk, eps_yield_init, sig_yield_init, modulus[0], s) + last_stress1;
				points[i].push_back(Point(xval, yval));
			}

			last_deform1 = points[i].back().x;
			last_stress1 = points[i].back().y;
		}


		write_to_file(points, FilePath + std::string("/points/model_transformed_points.alv"));

		std::cout << "Finished" << std::endl;
	}
	catch (const std::exception& e) { std::cout << e.what() << std::endl; }
}

void calc_init_modulus_distribution(const char* folder_path, vec& distribution) {
	distribution.resize(0);

	std::string folder = folder_path;
	std::string to_erase = folder + "\\";
	std::string check = ".alv";
	for (const auto & entry : std::filesystem::directory_iterator(folder)) {
		std::filesystem::path path = entry.path();
		std::string path_string = path.u8string();
		if (path_string.find(check) != std::string::npos) {
			/*size_t pos = path_string.find(to_erase);
			path_string.erase(pos, to_erase.length());*/

			//ищем модуль упругости
			dataPoint original_points;
			dataPoint transformed_points;
			vec modulus;
			dataDouble current_modulus;
			indexData pr_stress_ind;
			double average_modulus;

			read_data(original_points, transformed_points, path_string);
			calc_current_modulus(transformed_points, current_modulus);
			calc_proportional_stress(transformed_points, current_modulus, pr_stress_ind);
			calc_modulus(current_modulus, pr_stress_ind, modulus, average_modulus);

			distribution.push_back(modulus[0]);
		}
	}
}

void calc_init_yield_distribution(const char* folder_path, vecPoint& distribution) {
	distribution.resize(0);

	std::string folder = folder_path;
	std::string to_erase = folder + "\\";
	std::string check = ".alv";
	for (const auto & entry : std::filesystem::directory_iterator(folder)) {
		std::filesystem::path path = entry.path();
		std::string path_string = path.u8string();
		if (path_string.find(check) != std::string::npos) {
			/*size_t pos = path_string.find(to_erase);
			path_string.erase(pos, to_erase.length());*/

			//ищем модуль упругости
			dataPoint original_points;
			dataPoint transformed_points;
			vec modulus;
			dataDouble current_modulus;
			indexData pr_stress_ind;
			vec yield_stress, yield_deform;
			vec odquist;
			double average_modulus;

			read_data(original_points, transformed_points, path_string);
			calc_current_modulus(transformed_points, current_modulus);
			calc_proportional_stress(transformed_points, current_modulus, pr_stress_ind);
			calc_modulus(current_modulus, pr_stress_ind, modulus, average_modulus);

			calc_odqist(original_points, odquist);
			calc_yield_stress(transformed_points, modulus, pr_stress_ind, yield_stress, yield_deform);

			distribution.push_back(Point(yield_deform[0], yield_stress[0]));
		}
	}
}

void calc_odquist_given_modulus_distr(dataDouble& model_odquist, CopyUtility& s, vec& a, vec& b, vec& alpha, vec& odquist_medium, std::string distr_file_path) {

	//считаем из файла значени€
	vec modulus_data;
	std::ifstream inFile;
	inFile.open(distr_file_path);
	double x;
	while (inFile >> x)
		modulus_data.push_back(x);

	model_odquist.resize(modulus_data.size());

	inFile.close();

	/*CopyUtility s;
	vec a, b, w, odquist_medium;
	calc_everything(points_file_path, s, a, b, w, odquist_medium);*/

	for (int i = 0; i < modulus_data.size(); i++) {
		double E_init = modulus_data[i];

		calc_model_odqist_given_init_modulus(model_odquist[i], E_init, s, a, b, alpha, odquist_medium);

		std::cout << i + 1 << " sample of " << modulus_data.size() << std::endl;
	}
	std::cout << " Finished " << std::endl;
}

void calc_odquist_given_yield_distr(dataDouble& model_odquist, CopyUtility& s, vec& a, vec& b, vec& alpha, vec& odquist_medium, std::string deform_distr_file_path, const char* stress_distr_file_path) {

	//считаем из файла значени€
	vec yield_stress_data, yield_deform_data;
	std::ifstream inFilex, inFiley;
	inFiley.open(stress_distr_file_path);
	double y;

	while (inFiley >> y)
		yield_stress_data.push_back(y);
	
	inFiley.close();

	model_odquist.resize(yield_stress_data.size());

	/*CopyUtility s;
	vec a, b, w, odquist_medium;
	calc_everything(points_file_path, s, a, b, w, odquist_medium);*/

	for (int i = 0; i < yield_stress_data.size(); i++) {

		//найдем точку с таким напр€жением
		int curr_ind = 0;
		double curr_diff = INFINITY;
		for (int k = 1; k < s.transformed_points[0].size(); k++) {
			if (curr_diff > abs(yield_stress_data[i] - s.transformed_points[0][k].y)) {
				curr_ind = k;
				curr_diff = abs(yield_stress_data[i] - s.transformed_points[0][k].y);
			}
		}

		double yield_str_init;
		double yield_deform_init;

		//ближайша€ точка не последн€€
		bool flag1 = ((curr_ind + 1) != s.original_points[0].size());
		//точка между curr_ind и curr_ind + 1
		bool flag2 = ((s.original_points[0][curr_ind].y <= yield_stress_data[i]) && (yield_stress_data[i] <= s.original_points[0][curr_ind + 1].y));
		//ближайша€ точка не перва€
		bool flag3 = (curr_ind != 0);

		//когда ближайшие точки не перва€ и не последн€€
		if (flag1 && flag3) {
			if (flag2) {
				double x1 = s.transformed_points[0][curr_ind].y;
				double y1 = s.transformed_points[0][curr_ind].x;

				double x2 = s.transformed_points[0][curr_ind + 1].y;
				double y2 = s.transformed_points[0][curr_ind + 1].x;

				yield_str_init = yield_stress_data[i];
				yield_deform_init = (yield_str_init - x1) * (y2 - y1) / (x2 - x1) + y1;
			}
			else {
				double x1 = s.transformed_points[0][curr_ind].y;
				double y1 = s.transformed_points[0][curr_ind].x;

				double x2 = s.transformed_points[0][curr_ind - 1].y;
				double y2 = s.transformed_points[0][curr_ind - 1].x;

				yield_str_init = yield_stress_data[i];
				yield_deform_init = (yield_str_init - x1) * (y2 - y1) / (x2 - x1) + y1;
			}
		}
		else {
			//когда ближайша€ точка перва€
			if (curr_ind == 0) {
				double x1 = s.transformed_points[0][0].y;
				double y1 = s.transformed_points[0][0].x;

				double x2 = s.transformed_points[0][1].y;
				double y2 = s.transformed_points[0][1].x;

				yield_str_init = yield_stress_data[i];
				yield_deform_init = (yield_str_init - x1) * (y2 - y1) / (x2 - x1) + y1;
			}
			else {
				//ближайша€ последн€€ 
				//когда точка за пределами исходной кривой
				if (s.original_points[0].back().y < yield_stress_data[i]) {
					std::cout << "last yield point" << std::endl;
					yield_str_init = s.original_points[0][curr_ind].y;
					yield_deform_init = s.original_points[0][curr_ind].x;
				}
				else {
					//между последней и предпоследней точкой
					double x1 = s.transformed_points[0][curr_ind].y;
					double y1 = s.transformed_points[0][curr_ind].x;

					double x2 = s.transformed_points[0][curr_ind - 1].y;
					double y2 = s.transformed_points[0][curr_ind - 1].x;

					yield_str_init = yield_stress_data[i];
					yield_deform_init = (yield_str_init - x1) * (y2 - y1) / (x2 - x1) + y1;
				}

			}

		}		

		//double asa = s.transformed_points[0][s.yield_stress_ind[0]].x;
		//double asb = s.transformed_points[0][s.yield_stress_ind[0]].y;
		//calc_model_odqist_given_init_yield_stress(model_odquist[i], asa, asb, s, a, b, odquist_medium);
		calc_model_odqist_given_init_yield_stress(model_odquist[i], yield_deform_init, yield_str_init, s, a, b, alpha, odquist_medium);

		std::cout << i + 1 << " sample of " << yield_stress_data.size() << " " << yield_deform_init << " " << yield_str_init << std::endl;
	}
	std::cout << " Finished " << std::endl;
}

void calc_all(const char* name1) {
	try {
		_mkdir(name1);

		CopyUtility s1;
		vec a1, b1, alpha1, odquist_medium1;
		calc_everything(name1, s1, a1, b1, alpha1, odquist_medium1);

		//одквист по распределени€м
		std::string root1 = name1 + std::string("/odqist_data");
		_mkdir(root1.c_str());

		dataDouble model_odquist_E1, model_odquist_E1_with_crit;
		//indexData destruct_distro_E1;
		vec destruct_distro_E1;
		calc_odquist_given_modulus_distr(model_odquist_E1, s1, a1, b1, alpha1, odquist_medium1, "distro_data/random_modules.txt");
		//calc_with_destruction_crit(destruct_distro_E1, model_odquist_E1_with_crit, model_odquist_E1);
		calc_with_destruction_crit_continuous(destruct_distro_E1, model_odquist_E1_with_crit, model_odquist_E1);

		write_to_file(model_odquist_E1, name1 + std::string("/odqist_data/modelled_odquist_E.alv"));
		write_to_file(model_odquist_E1_with_crit, name1 + std::string("/odqist_data/modelled_odquist_E_wcrit.alv"));
		write_to_file(destruct_distro_E1, name1 + std::string("/odqist_data/E_wcrit_distro.alv"));


		dataDouble model_odquist_yield1, model_odquist_yield1_with_crit;
		//indexData destruct_distro_yield1;
		vec destruct_distro_yield1;
		calc_odquist_given_yield_distr(model_odquist_yield1, s1, a1, b1, alpha1, odquist_medium1, "distro_data/random_yield_deform.txt", "distro_data/random_yield_stress.txt");
		//calc_with_destruction_crit(destruct_distro_yield1, model_odquist_yield1_with_crit, model_odquist_yield1);
		calc_with_destruction_crit_continuous(destruct_distro_yield1, model_odquist_yield1_with_crit, model_odquist_yield1);

		write_to_file(model_odquist_yield1, name1 + std::string("/odqist_data/model_odqist_yield.alv"));
		write_to_file(model_odquist_yield1_with_crit, name1 + std::string("/odqist_data/model_odqist_yield_wcrit.alv"));
		write_to_file(destruct_distro_yield1, name1 + std::string("/odqist_data/yield_wcrit_distro.alv"));

	}
	catch (const std::exception& e) { std::cout << e.what() << std::endl; }
}

void calc_distros() {
	try {
		//распределени€
		vec modulus_distribution;
		calc_init_modulus_distribution("C:/Users/we/Documents/Diploma/Exp", modulus_distribution);
		write_to_file(modulus_distribution, "distro_data/modulus_distribution.alv");

		vecPoint yield_distribution;
		calc_init_yield_distribution("C:/Users/we/Documents/Diploma/Exp", yield_distribution);
		write_to_file_reverse(yield_distribution, "distro_data/yield_distribution.alv");
		write_to_file_second_coord(yield_distribution, "distro_data/yield_distribution_Y.alv");
	}
	catch (const std::exception& e) { std::cout << e.what() << std::endl; }
}

void merge_two_files(const char* file1_path, const char* file2_path, const char* destination) {

	dataPoint original_points;
	dataPoint transformed_points;

	read_data_no_scale(original_points, transformed_points, file1_path);
	int max_ind = original_points.size();

	std::ofstream outfile;

	outfile.open(destination/*, std::ios_base::app*/);
	for (int j = 0; j < original_points.size(); j++)
		for (int i = 0; i < original_points[j].size(); i++)
			outfile << j << "	" << original_points[j][i].y << "	" << original_points[j][i].x << std::endl;
	outfile.close();

	original_points.clear();
	transformed_points.clear();

	read_data_no_scale(original_points, transformed_points, file2_path);

	std::ofstream outfile1;

	outfile1.open(destination, std::ios_base::app);
	for (int j = 0; j < original_points.size(); j++)
		for (int i = 0; i < original_points[j].size(); i++)
			outfile1 << j + max_ind << "	" << original_points[j][i].y << "	" << original_points[j][i].x << std::endl;
	outfile1.close();

}

void calc_data_to_analyze() {
	std::vector<std::string> names;
	
	std::string name1 = "G39-22-1-7";
	std::string name2 = "G41-22-3-21";
	std::string name3 = "G42-22-1-17";
	std::string name4 = "G49-22-2-21";
	std::string name5 = "G66-22-8-1";
	std::string name6 = "G105-22-3-1";
	std::string name7 = "G112-22-1-13";
	std::string name8 = "G118-22-21";
	std::string name9 = "G131-22-6-7";
	std::string name10 = "G206-22-6-18";
	std::string name11 = "G219-22-6-9";
	std::string name12 = "G220-22-6-16";
	std::string name13 = "G223-22-6-13";
	std::string name14 = "G228-22-6-8";

	names.push_back(name1);
	names.push_back(name2);
	names.push_back(name3);
	names.push_back(name4);
	names.push_back(name6);
	names.push_back(name7);
	names.push_back(name8);
	names.push_back(name9);
	names.push_back(name10);
	names.push_back(name11);
	names.push_back(name12);
	names.push_back(name13);
	names.push_back(name14);

	std::string root1 =  std::string("to_analyze3/");
	_mkdir(root1.c_str());

#pragma omp parallel for num_threads(8)
	for (int i = 0; i < names.size(); i++) {
		try {

			const char* name = names[i].c_str();

			//_mkdir(name);

			CopyUtility s1;
			vec a1, b1, alpha1, odquist_medium1;
			calc_everything(name, s1, a1, b1, alpha1, odquist_medium1);

			double range = deformations_range(s1.original_points);
			std::ofstream outfile1;
			outfile1.open(root1 + std::to_string(i) + std::string("def_range.alv"));
			outfile1 << range << std::endl;
			outfile1.close();

			//одквист по распределени€м
			
			dataDouble model_odquist_E1, model_odquist_E1_with_crit;
			vec destruct_distro_E1;
			calc_odquist_given_modulus_distr(model_odquist_E1, s1, a1, b1, alpha1, odquist_medium1, "distro_data/random_modules.txt");
			calc_with_destruction_crit_continuous(destruct_distro_E1, model_odquist_E1_with_crit, model_odquist_E1);

			write_to_file(destruct_distro_E1,  root1 + std::to_string(i) + std::string("E_wcrit_distro.alv"));


			dataDouble model_odquist_yield1, model_odquist_yield1_with_crit;
			vec destruct_distro_yield1;
			calc_odquist_given_yield_distr(model_odquist_yield1, s1, a1, b1, alpha1, odquist_medium1, "distro_data/random_yield_deform.txt", "distro_data/random_yield_stress.txt");
			calc_with_destruction_crit_continuous(destruct_distro_yield1, model_odquist_yield1_with_crit, model_odquist_yield1);

			write_to_file(destruct_distro_yield1, root1 + std::to_string(i) + std::string("yield_wcrit_distro.alv"));

		}
		catch (const std::exception& e) { std::cout << e.what() << std::endl; }
	}
}

void calc_criteria_points() {

	std::vector<std::string> names;
	vecPoint points;

	std::string name1 = "G39-22-1-7";
	std::string name2 = "G41-22-3-21";
	std::string name3 = "G42-22-1-17";
	std::string name4 = "G49-22-2-21";
	std::string name5 = "G66-22-8-1";
	std::string name6 = "G105-22-3-1";
	std::string name7 = "G112-22-1-13";
	std::string name8 = "G118-22-21";
	std::string name9 = "G131-22-6-7";
	std::string name10 = "G206-22-6-18";
	std::string name11 = "G219-22-6-9";
	std::string name12 = "G220-22-6-16";
	std::string name13 = "G223-22-6-13";
	std::string name14 = "G228-22-6-8";

	names.push_back(name1);
	names.push_back(name2);
	names.push_back(name3);
	names.push_back(name4);
	names.push_back(name6);
	names.push_back(name7);
	names.push_back(name8);
	names.push_back(name9);
	names.push_back(name10);
	names.push_back(name11);
	names.push_back(name12);
	names.push_back(name13);
	names.push_back(name14);


	CopyUtility s;
	vec a, b, alpha, odquist_medium;
	calc_everything(names[0].c_str(), s, a, b, alpha, odquist_medium);
	write_to_file(s.odquist, "odquist2.alv");

//#pragma omp parallel for num_threads(8)
//	for (int i = 0; i < names.size(); i++) {
//		CopyUtility s;
//		vec a, b, alpha, odquist_medium;
//		calc_everything(names[i].c_str(), s, a, b, alpha, odquist_medium);
//
//		if(i == 3)
//			write_to_file(s.odquist, "odquist2.alv");
//
//		//double E_init = s.modulus[0];
//
//		//int additional_iterations = 0;
//
//		//model_odquist.resize(s.transformed_points.size() - 1 + additional_iterations);
//
//		////считаем параметр ќдквиста на каждом полуцикле (кроме последнего)
//		//vec odquist_gain;
//		//odquist_gain.resize(model_odquist.size());
//
//		//vec zero_def_elastic_x;
//		//zero_def_elastic_x.push_back(0.0);
//
//		////int yield_index_init = s.yield_stress_ind[0];
//		//double eps_yield_init = s.yield_deform[0];
//		//double sig_yield_init = s.yield_stress[0];
//
//		////значение деформации и напр€дежени€ в последней точке прошлого полуцикла
//		//double last_deform = s.transformed_points[0].back().x;
//		//double last_stress = s.transformed_points[0].back().y;;
//
//		////на первом полуцикле
//		//model_odquist[0] = 0.0;//s.odquist[0];
//		//odquist_gain[0] = 0.0;//-last_stress / (-E_init * interpolate_d(s.odquist[0], s.d, s.odquist, 1));//s.odquist[1] - s.odquist[0];
//
//
//
//		//double average_deform = 0.0;
//		//for (int i = 1; i < s.transformed_points.size(); i++)
//		//	average_deform += s.transformed_points[i].back().x;
//		//average_deform /= s.transformed_points.size();
//
//
//		//for (int i = 1; i < s.transformed_points.size() + additional_iterations - 1; i++) {
//		//	double accuracy = 1e-4;
//
//		//	//одквист
//		//	double eps1, eps2, gain;
//
//		//	double sign;
//		//	if (i % 2 != 0)
//		//		sign = -1.0;
//		//	else
//		//		sign = 1.0;
//
//		//	model_odquist[i] = model_odquist[i - 1] + odquist_gain[i - 1];
//
//		//	double odquistk = model_odquist[i];
//		//	double ak = interpolate_ab(odquistk, a, odquist_medium);
//		//	double bk = interpolate_ab(odquistk, b, odquist_medium);
//		//	double dk = interpolate_d(odquistk, s.d, s.odquist, i);
//		//	double final_deform;
//		//	/*if (i < s.transformed_points.size())
//		//		final_deform = s.transformed_points[i].back().x;
//		//	else
//		//		final_deform = average_deform;*/
//		//	final_deform = s.transformed_points[i].back().x;
//
//		//	//найдем точку пересечени€ осью абцисс продолжени€ упругих участков
//		//	//last_stress + sign * E_init * dk * deform = 0;
//		//	zero_def_elastic_x.push_back(-last_stress / (sign * E_init * dk));
//
//		//	odquist_gain[i] = abs(zero_def_elastic_x[i - 1] - zero_def_elastic_x[i]);
//
//
//		//	last_deform += sign * final_deform;
//		//	last_stress += sign * model_stress(final_deform, ak, bk, dk, eps_yield_init, sig_yield_init, E_init, s);
//
//		//}
//
//		points.push_back(Point(s.transformed_points.size() - 1, s.odquist.back()));
//	}

	//write_to_file(points, "criteria_points.alv");
}