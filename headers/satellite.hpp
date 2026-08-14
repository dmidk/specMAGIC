#pragma once

#include <stdexcept>
#include <string>
#include <array>
#include <cmath>
#include "constants.hpp"
#include "types.hpp"

namespace Satellite {

    MAGIC_EXACT calcObsTime(int start_hour, int start_minute, int line, int n_lines);

    int flipVertical(int line, int height);

    void geo2MTGImage(MAGIC_EXACT lat_rad, MAGIC_EXACT lon_rad,
        int nav_res, int col_off, int line_off,
        int max_cols, int max_lines, int& col, int& line);

}