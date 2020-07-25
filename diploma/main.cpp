#include "mainfunc.h"
#include "statistics.h"

//#include <windows.h>

int main() {
	setlocale(LC_ALL, "Russian");

	//calc_criteria_points();

	//char name1[] = "G219-22-6-9";
	//char name2[] = "G105-22-3-1";
	//char name3[] = "G131-22-6-7";
	//char name4[] = "G42-22-1-17";

	//calc_distros();

	//calc_all(name1);
	//calc_all(name2);
	//calc_all(name3);
	//calc_all(name4);

	//merge_two_files("exp_second/G47-22-2-21.alv", "exp_second/G49-22-2-21_continue.alv", "exp_second/G49-22-2-21MERGED.alv");


	calc_data_to_analyze();


	/*std::vector<std::string> names;

	std::string name1 = "G39-22-1-7";
	std::string name2 = "G41-22-3-21";
	std::string name3 = "G42-22-1-17";
	std::string name4 = "G105-22-3-1";
	std::string name5 = "G112-22-1-13";
	std::string name6 = "G118-22-21";
	std::string name7 = "G131-22-6-7";
	std::string name8 = "G206-22-6-18";
	std::string name9 = "G219-22-6-9";
	std::string name10 = "G220-22-6-16";
	std::string name11 = "G223-22-6-13";
	std::string name12 = "G228-22-6-8";
	std::string name13 = "G49-22-2-21";
	std::string name14 = "G66-22-8-1";

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


	for (int i = 0; i < names.size(); i++) {
		CopyUtility s1;
		vec a1, b1, odquist_medium1, alpha1;
		calc_everything(names[i].c_str(), s1, a1, b1, alpha1, odquist_medium1);
	}*/
	
	/*CopyUtility s1;
	vec a1, b1, odquist_medium1, alpha1;
	calc_everything("G219-22-6-9", s1, a1, b1, alpha1, odquist_medium1);*/


	//test_all(12, "to_analyze3");

	getchar();
	return 0;
}


