#pragma once
#include "halfcycle.h"
#include "fatigue.h"
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <algorithm>
#include <direct.h>
#include "func.h"
#include <omp.h>

void calc_everything(const char* FilePath, CopyUtility& s, vec& a, vec& b, vec& alpha, vec& odquist_medium);

void calc_init_modulus_distribution(const char* folder_path, vec& distribution);
void calc_init_yield_distribution(const char* folder_path, vecPoint& distribution);

void calc_odquist_given_modulus_distr(dataDouble& model_odquist, CopyUtility& s, vec& a, vec& b, vec& alpha, vec& odquist_medium, std::string distr_file_path);
void calc_odquist_given_yield_distr(dataDouble& model_odquist, CopyUtility& s, vec& a, vec& b, vec& alpha, vec& odquist_medium, std::string deform_distr_file_path, const char* stress_distr_file_path);

void calc_distros();

void calc_all(const char* FilePath);

void merge_two_files(const char* file1_path, const char* file2_path, const char* destination);

void calc_data_to_analyze();

void calc_criteria_points();



