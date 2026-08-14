#pragma once 
#include <netcdf.h>
#include <string>
#include <cstring> 
#include <vector>
#include <cstdio>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include "netcdf.hpp"
#include "navigation.hpp"
#include "matrix.hpp"

std::string makeFilename(DateTime& timestamp,  std::string& variable_name, std::string& out_path);

int writeToNetcdf(const std::string& filename, Geography& geo, Matrix& box, int year, int month, int day,
    double seconds_since_midnight, std::string name);