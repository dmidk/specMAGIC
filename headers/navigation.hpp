#pragma once 
#include <array>
#include <algorithm>
#include <fstream>
#include <cassert>
#include <cctype>
#include <cstdio>  
#include <utility>  // std::move
#include <math.h> // floor
#include <cmath>
#include <stdexcept>

#include "constants.hpp"
#include "helpers.hpp"
#include "types.hpp"

struct DateTime {
    
    unsigned int year   = 0;
    unsigned int month  = 0;
    unsigned int day    = 0;
    unsigned int hour   = 0;
    unsigned int minute = 0;

    // Parameterised ctor
    DateTime(int yr, int mth, int dy, int hr, int min) : year(yr), month(mth), 
        day(dy), hour(hr), minute(min) {}

    // Define default ctor, otherwise build errors...
    DateTime() : year(0), month(0), day(0), hour(0), minute(0) {}

    // Idiot check
    bool valid() {

        bool is_valid = false;
        if (year > 2000 && month <= 12 && day <= 31 && hour <= 23 && minute <= 59) is_valid = true;

        return is_valid;
    }
};

// Teeny tiny struct just to contain the variables
// So they're not all wobbling around
// Climatology etc. params for a given pixel
struct PixelClimate {

    MAGIC_REAL aod;
    MAGIC_REAL ssa;
    MAGIC_REAL gg;
    MAGIC_REAL h2o; 
    MAGIC_REAL o3;
    MAGIC_EXACT surface_albedo;
    MAGIC_EXACT albedo_correction;

    // Parameterised ctor for convenience
    PixelClimate(MAGIC_REAL ao, MAGIC_REAL ss, MAGIC_REAL g, MAGIC_REAL h2, MAGIC_REAL o, MAGIC_EXACT alb, MAGIC_EXACT alb_corr) : aod(ao), ssa(ss), 
        gg(g), h2o(h2), o3(o), surface_albedo(alb), albedo_correction(alb_corr) {}

    // Required to define default ctor otherwise build errors
    PixelClimate() : aod(0), ssa(0), gg(0), h2o(0), o3(0), surface_albedo(0), albedo_correction(0) {}


};


struct Geography {

    // Number of points in each dimension
    unsigned int nlat = 0;
    unsigned int nlon = 0;

    // Arrays of actual coordinate points
    // in RADIANS
    MAGIC_EXACT *lat = nullptr;
    MAGIC_EXACT *lon = nullptr;

    // NO default ctor
    // Any instance should be immediately initialised
    // Use the parameterised ctor

    // Constructor makes the allocations
    Geography(unsigned int nlat_in, unsigned int nlon_in)
        : nlat(nlat_in), nlon(nlon_in) {
        lat = new MAGIC_EXACT[nlat];
        lon = new MAGIC_EXACT[nlon];
    }

    // Destructor: because we use new
    ~Geography() {
        delete[] lat;
        delete[] lon;
    }

    // No copy because I hate it
    Geography(const Geography&) = delete;
    Geography& operator=(const Geography&) = delete;

    // Allow moves
    Geography(Geography&& other) noexcept {
        *this = std::move(other);
    }

    // Move assignment
    Geography& operator=(Geography&& other) noexcept {
        if (this != &other) {
            delete[] lat;
            delete[] lon;

            nlat = other.nlat;
            nlon = other.nlon;
            lat  = other.lat;
            lon  = other.lon;

            other.nlat = other.nlon = 0;
            other.lat = other.lon = nullptr;
        }
        return *this;
    }

    void makeCoordinates(Config c) {

        MAGIC_REAL latbegin = c.latbegin;
        MAGIC_REAL lonbegin = c.lonbegin;
        MAGIC_REAL dxy = c.dxy;
        MAGIC_REAL deltalon = c.deltalon;

        // Make latitude coords
        for (unsigned int l = 0; l < nlat; l++) {
            lat[l] = static_cast<MAGIC_EXACT>(latbegin + dxy * l);
        }

        // Make longitude coords
        for (unsigned int l = 0; l < nlon; l++) {
            lon[l] = static_cast<MAGIC_EXACT>(lonbegin + dxy * l + deltalon);

            // We want -180 to 180, not 0 to 360
            if (lon[l] > 180.) {
                lon[l] -= 360.;
            }
        }
    }
};



struct Area {

    // This is NOT THE SAME as Geography struct above
    // Geography is global, for the entire image
    // This is for the loop iterations, i.e. per pixel

    // These are the grid indices of the image
    unsigned int nlat;   
    unsigned int nlon; 

    // Geographical lat/lon, in degrees
    MAGIC_EXACT degrees_lat;
    MAGIC_EXACT degrees_lon;

    // Geographical lat/lon, in radians
    MAGIC_EXACT lat;
    MAGIC_EXACT lon;
    MAGIC_EXACT deltalon_rad; 

    Area() : nlat(0), nlon(0), degrees_lat(0), degrees_lon(0), lat(0), lon(0) {}

    Area(unsigned nlat_in, unsigned nlon_in, MAGIC_EXACT deglat, MAGIC_EXACT deglon, MAGIC_EXACT lat_in, MAGIC_EXACT lon_in, MAGIC_EXACT deltalon_rad_in) : 
        nlat(nlat_in), nlon(nlon_in), degrees_lat(deglat), degrees_lon(deglon), 
        lat(lat_in), lon(lon_in), deltalon_rad(deltalon_rad_in) {}

    void makeArea(const Geography& c, int iter, MAGIC_EXACT delta_lon) {

        nlat = iter / c.nlon;
        nlon = iter % c.nlon;

        degrees_lat = c.lat[nlat];
        degrees_lon = c.lon[nlon];

        lat = degrees_lat / RADTODEGREE;
        lon = degrees_lon / RADTODEGREE;
        deltalon_rad = delta_lon / RADTODEGREE;

        assert(nlat >= 0 && nlat < c.nlat);
        assert(nlon >= 0 && nlon < c.nlon);
    }
};


