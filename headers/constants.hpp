#pragma once 
#include "types.hpp"

constexpr int NAV_CORRECTION_COLOFFSET = 1;      // empirically determined correction on column offset for geolocation problems 
constexpr int NAV_CORRECTION_LINOFFSET = 0;      // empirically determined correction on line offset for geolocation problems 

// Radiation bands
// Indices of first and last kato bands
constexpr int KATO_MIN = 1; 
constexpr int KATO_MAX = 32;

// Useful maths
constexpr MAGIC_EXACT RADTODEGREE = 57.29577951308232087540;
constexpr MAGIC_EXACT PI = 3.141592653589793;

/// Satellite distance from Earth center [km]
constexpr MAGIC_EXACT SATELLITE_RADIUS_KM = 42164.0;

/// Earth polar radius [km]
constexpr MAGIC_EXACT EARTH_POLAR_RADIUS_KM = 6356.5358;

/// Ellipsoid parameters
constexpr MAGIC_EXACT EPSI2 = 0.00676765;     // (Re^2 - Rp^2) / Re^2
constexpr MAGIC_EXACT RPE2  = 0.99323235;     // (Rp / Re)^2

// For computing the line offset per pixel
constexpr MAGIC_EXACT FULL_DISK_SCAN_MIN = 9.0;  // is this accurate for dk? TODO: Check
constexpr MAGIC_EXACT HOURS_PER_MIN = 1.0 / 60.0;
constexpr MAGIC_EXACT COS_SOLAR_ZENITH_MIN = 1e-3;

// For sun geometry
constexpr MAGIC_EXACT HOURS_AT_NOON = 12.0;
constexpr MAGIC_EXACT DEG_PER_HOUR  = 15.0;       // Earth rotation
constexpr MAGIC_EXACT RAD_PER_DEG   = PI / 180.0;
constexpr MAGIC_EXACT MIN_PER_HOUR  = 60.0;

// About 67 degrees
// Close to limb
// Limit of what we can reliably use
constexpr MAGIC_EXACT MAX_VALID_LIMB_RADIUS = 1.2;     // in radians

// Cloud index thresholds
constexpr MAGIC_EXACT CAL_MIN_CLEAR = -0.05;     // below this: very clear / dark
constexpr MAGIC_EXACT CAL_NOISE_THRESHOLD = 0.04; // below this: treat as clear
constexpr MAGIC_EXACT CAL_LINEAR_MAX = 0.9;      // end of linear regime
constexpr MAGIC_EXACT CAL_TRANSITION_MAX = 1.03; // end of transition regime

// Cloud transmission coefficients
constexpr MAGIC_EXACT K_CLEAR_EXTREME = 1.05;
constexpr MAGIC_EXACT K_OVERCAST = 0.05;

// Clear sky index coefficients
constexpr MAGIC_EXACT K_A0 = 1.115;
constexpr MAGIC_EXACT K_A1 = -1.781;
constexpr MAGIC_EXACT K_A2 = 0.726;

// Helpful for val scaling
constexpr MAGIC_EXACT CAL_SCALE = 100.0;
constexpr int CAL_MIN_INT = 0;
constexpr int CAL_MAX_INT = 100;
constexpr int CAL_UNDEFINED = -1;

// For the CAL calculation
constexpr MAGIC_EXACT DARK_OFFSET = 0;




