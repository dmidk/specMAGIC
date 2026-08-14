#pragma once
#include "matrix.hpp"
#include "write.hpp"
#include "types.hpp"
#include "navigation.hpp"

// Scratch memory for spectral arrays
// MUST be one per thread
struct Scratch {

    int n = 0;        // number of elements

    // Whether we are having a valid daytime pixel
    bool valid = false;

    // The same as out_G and out_B in DWD code
    MAGIC_REAL* DNI_spectral = nullptr;
    MAGIC_REAL* GHI_spectral = nullptr;

    MAGIC_REAL kB = 0.;
    MAGIC_REAL kG = 0.;


    // Constructor ALWAYS allocs out arrays
    Scratch(int n_in)
        : n(n_in) {
        DNI_spectral = new  MAGIC_REAL[n];
        GHI_spectral = new  MAGIC_REAL[n];

        for (int i = 0; i < n; i++) {
            DNI_spectral[i] = 0.;
            GHI_spectral[i] = 0.;
        }
    }

    // Destructor, because we use new
    ~Scratch() {
        delete[] DNI_spectral;
        delete[] GHI_spectral;

    }

    // No copy
    Scratch(const Scratch&) = delete;
    Scratch& operator=(const Scratch&) = delete;

    // Allow moves
    Scratch(Scratch&& other) noexcept {
        *this = std::move(other);
    }

    Scratch& operator=(Scratch&& other) noexcept {
        if (this != &other) {
            delete[] DNI_spectral;
            delete[] GHI_spectral;

            n  = other.n;
            DNI_spectral  = other.DNI_spectral;
            GHI_spectral  = other.GHI_spectral;
            kB = other.kB;
            kG = other.kG;

            other.n = 0;
            other.DNI_spectral = nullptr;
            other.GHI_spectral = nullptr;
        }
        return *this;
    }

    void reset() {
     for (int i = 0; i < n; i++) {
            DNI_spectral[i] = 0.;
            GHI_spectral[i] = 0.;
        }
        
    }
};

// The big radiation matrices to store all data
struct Output {
    Matrix GHI;     // Surface incoming shortwave
    Matrix DNI;     // Direct normalised irradiance
    Matrix CAL;     // Effective cloud albedo
    Matrix CSR;     // Clear-sky Radiance

    Output(int rows, int cols)
        : GHI(rows, cols),
          DNI(rows, cols),
          CAL(rows, cols), 
          CSR(rows, cols) {}

    void toFile(DateTime& timestamp, Geography& geo, Config& c) {

        double seconds = timestamp.hour * 3600.0 +
            timestamp.minute * 60.0;

        std::string variable_name = "GHI";
        std::string out = c.path+c.out_global_path;
        std::string fname = makeFilename(timestamp, variable_name, out);

        writeToNetcdf(fname, geo, GHI, timestamp.year, timestamp.month, timestamp.day,
            seconds, variable_name);

        variable_name = "DNI";
        out = c.path+c.out_beam_path;
        fname = makeFilename(timestamp, variable_name, out);
        writeToNetcdf(fname, geo, DNI, timestamp.year, timestamp.month, timestamp.day,
            seconds, variable_name);

        variable_name = "CAL";
        out = c.path+c.out_cal_path;
        fname = makeFilename(timestamp, variable_name, out);
        writeToNetcdf(fname, geo, CAL, timestamp.year, timestamp.month, timestamp.day,
            seconds, variable_name);

        variable_name = "CSR";
        out = c.path+c.out_clear_path;
        fname = makeFilename(timestamp, variable_name, out);
        writeToNetcdf(fname, geo, CSR, timestamp.year, timestamp.month, timestamp.day,
            seconds, variable_name);



    }

    // Could do overloads here instead, but not much point

    inline void encodeCAL(MAGIC_REAL cal, Area a) {

        // Scale to output format
        MAGIC_INT scaled = static_cast<MAGIC_INT>(std::rint(CAL_SCALE * cal));

        // Clamp to valid range
        if (scaled < CAL_MIN_INT) scaled = CAL_MIN_INT;
        if (scaled > CAL_MAX_INT) scaled = CAL_MAX_INT;

        this->CAL(a.nlat, a.nlon) = static_cast<MAGIC_INT>(scaled);
    }

    inline void encodeDNI(MAGIC_REAL dni, MAGIC_REAL angle, Area a) {

        // Div by solar zenith angle
        if (dni >= 0 && angle >= 0) dni /= angle;
        
        this->DNI(a.nlat, a.nlon) = static_cast<MAGIC_INT>(std::lround(dni));
    }


    inline void encodeGHI(MAGIC_REAL ghi, Area a) {
        this->GHI(a.nlat, a.nlon) = static_cast<MAGIC_INT>(std::lround(ghi));

    }

    inline void encodeCSR(MAGIC_REAL csr, Area a) {
        this->CSR(a.nlat, a.nlon) = static_cast<MAGIC_INT>(std::lround(csr));
    }

};