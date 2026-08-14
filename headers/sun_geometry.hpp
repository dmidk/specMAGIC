#pragma once 
#include <array>
#include "constants.hpp"
#include "navigation.hpp"
#include "types.hpp"

// Tiny helper for all the sun geom info
struct SolarParameters {
    MAGIC_REAL gamma;        // fractional year angle [rad]
    MAGIC_REAL distance;     // Earth–Sun distance factor
    MAGIC_REAL declination;  // solar declination [rad]
    MAGIC_REAL eq_time_min;  // equation of time [minutes]
    MAGIC_REAL sza;         // solar zenith angle [rad]
    MAGIC_REAL cos_sza;     // cosine of the solar zenith angle
    MAGIC_REAL subsolar_lon;    // subsolar longitude [radians]
    MAGIC_REAL corrected_cos_sza;   // with corrections if large SZA

};

namespace SunGeometry {

    PixelState calcPixelValidity(MAGIC_EXACT nadir, MAGIC_REAL cos_sza);

    MAGIC_REAL subsolarLongitude(MAGIC_EXACT t_utc_hours, MAGIC_EXACT eq_time_min);

    int dayOfYear(int day, int month, int year);

    MAGIC_REAL cosSolarZenith (MAGIC_EXACT lat_rad, MAGIC_EXACT lon_rad,
        MAGIC_REAL declination_rad, MAGIC_REAL subsolar_lon_rad);

    SolarParameters solarParameters(DateTime timestamp, MAGIC_EXACT UTC_time, Area a);
    
    void correctZenithAngle(SolarParameters& s);
    

}