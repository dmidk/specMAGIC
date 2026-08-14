#include "../headers/sun_geometry.hpp"

namespace SunGeometry {

    PixelState calcPixelValidity(MAGIC_EXACT nadir, MAGIC_REAL cos_sza) {
        
        PixelState state;

        if (nadir > MAX_VALID_LIMB_RADIUS) {
            state = PixelState::Invalid;
        }
        else if (cos_sza > (MAGIC_REAL) 0.0 && cos_sza < (MAGIC_REAL) COS_SOLAR_ZENITH_MIN) {
            state = PixelState::Invalid;   // twilight or otherwise problematic
        }
        else if (cos_sza >= (MAGIC_REAL) COS_SOLAR_ZENITH_MIN) {
            state = PixelState::Day;
        }
        else {
            state = PixelState::Night;
        }

        return state;


    }

    // Does what it says on the tin
    MAGIC_REAL subsolarLongitude(MAGIC_EXACT t_utc_hours, MAGIC_EXACT eq_time_min) {

        // True solar time at Greenwich [hours]
        const MAGIC_EXACT true_solar_time =
            t_utc_hours + eq_time_min / MIN_PER_HOUR;

        // Longitude where Sun is overhead
        // in RADIANS
        return static_cast<MAGIC_REAL>( (HOURS_AT_NOON - true_solar_time)
            * DEG_PER_HOUR
            * RAD_PER_DEG) ;
    }


    /// Compute day of year (1–365 or 366)
    int dayOfYear(int day, int month, int year) {
        
        static constexpr std::array<int, 12> lastDayOfPrevMonth{
            0, 31, 59, 90, 120, 151,
            181, 212, 243, 273, 304, 334
        };

        int doy = lastDayOfPrevMonth[month - 1] + day;

        // Account for leap years
        if ((year % 4 == 0) && month > 2) {
            doy += 1;
        }

        return doy;
    }

    // Compute cosine of solar zenith angle
    MAGIC_REAL cosSolarZenith (MAGIC_EXACT lat_rad, MAGIC_EXACT lon_rad,
            MAGIC_REAL declination_rad, MAGIC_REAL subsolar_lon_rad) {
        
        MAGIC_REAL hour_angle = (MAGIC_REAL) lon_rad - subsolar_lon_rad;

        MAGIC_REAL cos_theta =
            std::sin((MAGIC_REAL) lat_rad) * std::sin(declination_rad) +
            std::cos((MAGIC_REAL) lat_rad) * std::cos(declination_rad) * std::cos(hour_angle);

        return cos_theta;

    }

    /// Get solar parameters from day of year
    SolarParameters solarParameters(DateTime timestamp, MAGIC_EXACT UTC_time, Area a) {

        // get the day of year
        int doy = SunGeometry::dayOfYear(timestamp.day, timestamp.month, timestamp.year);

        SolarParameters s{};

        // Fractional year angle [rad]
        s.gamma = 2.0 * PI * (doy - 1.0) / 365.0;

        const MAGIC_REAL cg  = std::cos(s.gamma);
        const MAGIC_REAL sg  = std::sin(s.gamma);
        const MAGIC_REAL c2g = std::cos(2.0 * s.gamma);
        const MAGIC_REAL s2g = std::sin(2.0 * s.gamma);

        // Earth–Sun distance correction factor
        s.distance = static_cast<MAGIC_REAL>(
            1.00010 +
            0.034221 * cg +
            0.001280 * sg +
            0.000719 * c2g +
            0.000077 * s2g ) ;

        // Solar declination [rad]
        s.declination = static_cast<MAGIC_REAL>(
            0.006918 -
            0.399912 * cg +
            0.070257 * sg -
            0.006758 * c2g +
            0.000907 * s2g -
            0.002697 * std::cos(3.0 * s.gamma) +
            0.001480 * std::sin(3.0 * s.gamma) );

        // Equation of time [minutes]
        s.eq_time_min = static_cast<MAGIC_REAL>( 
            (0.000075 +
            0.001868 * cg -
            0.032077 * sg -
            0.014615 * c2g -
            0.040890 * s2g)
            * (4.0 * 180.0 / PI) );

        s.subsolar_lon = subsolarLongitude(UTC_time, s.eq_time_min);

        s.cos_sza = cosSolarZenith(a.lat, a.lon,
                        s.declination, s.subsolar_lon);

        s.sza = acos(s.cos_sza);

        return s;
    }

    // If the pixel has a very large SZA, this can cause RMIN to blow up
    // This is a correction for an FCI measuring fault
    void correctZenithAngle(SolarParameters& s) {

        // FCI correction for measuring fault 

        s.corrected_cos_sza = 0; 

        if (s.sza > (MAGIC_REAL) 1.4) s.corrected_cos_sza = s.cos_sza + (MAGIC_REAL) 0.6 * (s.sza - (MAGIC_REAL) 1.4);
        else s.corrected_cos_sza = s.cos_sza;

    }



}   // Namespace SunGeometry
