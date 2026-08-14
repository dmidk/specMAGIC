#include "../headers/satellite.hpp"

namespace Satellite {

    // Gets the exact obs time for a given line 
    // As the satellite takes some time to take the entire picture
    // TODO: The old code may have been wrapping the time incorrectly. Check this. 
    // TODO: Get this directly from the filename?
    MAGIC_EXACT calcObsTime(int start_hour, int start_minute, int line, int n_lines) {
        
        MAGIC_EXACT observation_time_utc =
            start_hour +
            start_minute / 60.0 +
            MAGIC_EXACT(line) * (FULL_DISK_SCAN_MIN / n_lines / 60.0);
            
        if (observation_time_utc < 0.0) return 24.0 - observation_time_utc;

        return observation_time_utc;

        // This seems more correct but does not produce same results as specmagic now
        // return fmod(observation_time_utc + 24.0, 24.0);

    }

    // Depending on how image origin is defined this may be needed
    int flipVertical(int line, int height) {
        return height - line;
    }


    /**
     * @brief Convert geographic coordinates (lat, lon) to MTG image pixel coordinates.
     *
     * This function implements the official geostationary projection used by MTG
     * (Meteosat Third Generation). It maps a geodetic Earth position (latitude,
     * longitude) to line/column indices in the MTG image.
     *
     * The calculation:
     *  - uses a WGS84-like Earth ellipsoid
     *  - accounts for satellite height and Earth curvature
     *  - checks satellite visibility (Earth limb)
     *
     * @param lat_rad   Latitude in radians
     * @param lon_rad   Longitude in radians (relative to satellite sub-point)
     * @param nav_res   Navigation resolution (e.g. 222 or 667 for MTG products)
     * @param col_off   Column offset (image center)
     * @param line_off  Line offset (image center)
     * @param max_cols  Number of columns in the image
     * @param max_lines Number of lines in the image
     * @param[out] col  Resulting column index (0-based)
     * @param[out] line Resulting line index (0-based)
     *
     * @throws std::runtime_error if the point is not visible or outside the image
     */
    void geo2MTGImage(MAGIC_EXACT lat_rad, MAGIC_EXACT lon_rad,
        int nav_res, int col_off, int line_off,
        int max_cols, int max_lines,
        int& col, int& line) {

        // -----------------------------
        // MTG navigation resolution
        // -----------------------------

        // Not currently in use!

        //MAGIC_EXACT navigation_resolution;   // [urad per pixel]
        int resolution_factor;

        if (nav_res == 222) {
            // navigation_resolution = 222.623596;
            resolution_factor = 3;
        }
        else if (nav_res == 667) {
            // navigation_resolution = 667.2044067;
            resolution_factor = 1;
        }
        else {
            throw std::runtime_error(
                "Invalid MTG navigation resolution: " + std::to_string(nav_res));
        }

        // -------------------------------------------------
        // Convert geodetic latitude -> geocentric latitude
        // -------------------------------------------------

        const MAGIC_EXACT geocentric_lat =
            std::atan(RPE2 * std::tan(lat_rad));

        const MAGIC_EXACT cos_lat = std::cos(geocentric_lat);
        const MAGIC_EXACT sin_lat = std::sin(geocentric_lat);

        // --------------------------------------------
        // Earth radius at given latitude (ellipsoid)
        // --------------------------------------------

        const MAGIC_EXACT earth_radius =
            EARTH_POLAR_RADIUS_KM /
            std::sqrt(1.0 - EPSI2 * cos_lat * cos_lat);

        // --------------------------------------------
        // Vector from satellite to Earth surface point
        // --------------------------------------------

        const MAGIC_EXACT r1 =
            SATELLITE_RADIUS_KM -
            earth_radius * cos_lat * std::cos(lon_rad);

        const MAGIC_EXACT r2 =
        -earth_radius * cos_lat * std::sin(lon_rad);

        const MAGIC_EXACT r3 =
            earth_radius * sin_lat;

        const MAGIC_EXACT range =
            std::sqrt(r1 * r1 + r2 * r2 + r3 * r3);

        // --------------------------------------------
        // Viewing angles (satellite scan angles)
        // --------------------------------------------

        const MAGIC_EXACT alpha = std::atan2(r2, r1);   // east-west scan angle
        const MAGIC_EXACT beta  = std::asin(r3 / range); // north-south scan angle

        // --------------------------------------------
        // Visibility check (Earth limb test)
        // --------------------------------------------

        const MAGIC_EXACT visibility =
            SATELLITE_RADIUS_KM * std::cos(alpha) * std::cos(beta) -
            range * (std::pow(std::cos(beta), 2.0) +
                    RPE2 * std::pow(std::sin(beta), 2.0));

        if (visibility < 0.0) {
            throw std::runtime_error(
                "Geographic point is not visible from MTG satellite");
        }

        // --------------------------------------------
        // Convert angles to image coordinates
        // --------------------------------------------

        constexpr MAGIC_EXACT MTG_SCAN_SCALE = 3712.0 / 17.832;

        const MAGIC_EXACT column =
            col_off -
            alpha * resolution_factor * MTG_SCAN_SCALE * (180.0 / PI);

        const MAGIC_EXACT line_val =
            line_off -
            beta * resolution_factor * MTG_SCAN_SCALE * (180.0 / PI);

        // --------------------------------------------
        // Round to nearest pixel
        // --------------------------------------------

        col  = static_cast<int>(std::lround(column));
        line = static_cast<int>(std::lround(line_val));

        // --------------------------------------------
        // Image bounds check
        // --------------------------------------------

        if (col < 0 || col >= max_cols ||
            line < 0 || line >= max_lines) {
            throw std::runtime_error(
                "Mapped pixel lies outside MTG image bounds");
        }
    }

}   // Namespace Navigation

