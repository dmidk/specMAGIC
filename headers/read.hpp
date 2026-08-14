#pragma once
#include <string>
#include <type_traits> // std::underlying_type_t
#include <fstream>
#include <stdexcept>
#include <string>
#include <charconv>
#include <cctype>
#include "helpers.hpp"

using namespace std;

Config loadConfig(const string& path_to_config, const string& path_to_home);