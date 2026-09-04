#pragma once
#include <string>
#include <fstream>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cstddef>
#include <cmath>
#include <iostream>
#include "types.hpp"
#include "constants.hpp"
#include "navigation.hpp"


namespace Tables {

    struct GroundAlbedo {
        static constexpr int Lands = 20; 

        int bands;  
        MAGIC_REAL* data; 
        MAGIC_REAL** a; 

        GroundAlbedo(int bandsdim) : bands(bandsdim), 
            data(nullptr), a(nullptr) {
            
            data = new MAGIC_REAL [Lands * bands];

            a = new MAGIC_REAL*[Lands];
            for (int land = 0; land < Lands; land++) {
                a[land] = data + land * bands;
            }
            
        }

        ~GroundAlbedo() {
            if (a != nullptr) delete [] a;
            delete[] data;
        }

        // Disable copy constructors so we don't go insane later
        GroundAlbedo(const GroundAlbedo&) = delete;
        GroundAlbedo& operator=(const GroundAlbedo&) = delete;

        // File format: first line is header (ignored), then MAGIC_REALs.
        // Spectral : 20 * bands MAGIC_REALs total
        bool read(std::string& filepath) {
            std::ifstream in(filepath);
            if (!in) return false;

            // Skip header 
            std::string header;
            getline(in, header);

            // Read the data
            for (int land = 0; land < Lands; land++) {
                for (int band = 0; band < bands; band++) {
                    if (!(in >> data[land * bands + band])) {
                        return false; // parse error or not enough numbers
                    }
                }
            }
            return true;
        }

        void fromFile(Config c) {

            std::string filepath = c.path + c.gr_alb_file;

            bool success = read(filepath);
            if (!success) {
                printf("Error reading ground albedo file! \n");

            }

        }
        
        MAGIC_REAL getAlbedo(int land, int band) {
            if (land < 0 || land >= Lands) {
                printf("Invalid land index! %d\n", land);
                return 0.0;
            }

            return a[land][band];
        }


    };

    struct RTM {
        // Dimensions
        int ggdim;
        int ssadim;
        int aoddim;
        int bandsdim;

        // Axis arrays
        MAGIC_REAL* lambda;  // [bandsdim]
        MAGIC_REAL* aod;     // [aoddim]
        MAGIC_REAL* ssa;     // [ssadim]
        MAGIC_REAL* gg;      // [ggdim]

        // 4D fields: [ggdim][ssadim][aoddim][bandsdim]
        MAGIC_REAL**** Im;
        MAGIC_REAL**** gtau;
        MAGIC_REAL**** ag;
        MAGIC_REAL**** btau;
        MAGIC_REAL**** ab;

        // Internal storage (so we can delete correctly)
        // For each 4D field we store:
        // - base pointer layers (p4, p3, p2) and contiguous data.
        struct Field4D {
            MAGIC_REAL**** p4 = nullptr;
            MAGIC_REAL***  p3 = nullptr;
            MAGIC_REAL**   p2 = nullptr;
            MAGIC_REAL*    data = nullptr;
        };

        Field4D Im_mem, gtau_mem, ag_mem, btau_mem, ab_mem;

        RTM(int ggdim_, int ssadim_, int aoddim_, int bandsdim_)
            : ggdim(ggdim_), ssadim(ssadim_), aoddim(aoddim_), bandsdim(bandsdim_),
            lambda(nullptr), aod(nullptr), ssa(nullptr), gg(nullptr),
            Im(nullptr), gtau(nullptr), ag(nullptr), btau(nullptr), ab(nullptr) {
            // Axis allocations
            lambda = new MAGIC_REAL[bandsdim];
            aod    = new MAGIC_REAL[aoddim];
            ssa    = new MAGIC_REAL[ssadim];
            gg     = new MAGIC_REAL[ggdim];

            // 4D field allocations
            Im   = allocate4D(Im_mem);
            gtau = allocate4D(gtau_mem);
            ag   = allocate4D(ag_mem);
            btau = allocate4D(btau_mem);
            ab   = allocate4D(ab_mem);
        }

        ~RTM() {
            // Free fields
            free4D(Im_mem);
            free4D(gtau_mem);
            free4D(ag_mem);
            free4D(btau_mem);
            free4D(ab_mem);

            // Free axes
            delete[] gg;
            delete[] ssa;
            delete[] aod;
            delete[] lambda;
        }

        // No copying (raw owning pointers)
        RTM(const RTM&) = delete;
        RTM& operator=(const RTM&) = delete;

        // --- Read the aerosol LUT file ---
        // For SSA in [0..ssadim-1]:
        //   For GG in [0..ggdim-1]:
        //     header line like: "ssa= 0.70 gg= 0.60"
        //     For AOD in [0..aoddim-1]:
        //       For L in [0..bandsdim-1]:
        //         7 MAGIC_REALs: aod lambda Im gtau ag btau ab
        bool read(const std::string& filepath) {
            std::ifstream in(filepath);
            if (!in) return false;

            for (int SSA = 0; SSA < ssadim; SSA++) {
                for (int GG = 0; GG < ggdim; GG++) {

                    // Parse header tokens: "ssa=" value "gg=" value
                    std::string tok1, tok2;
                    if (!(in >> tok1 >> ssa[SSA] >> tok2 >> gg[GG])) {
                        return false;
                    }

                    for (int AOD = 0; AOD < aoddim; AOD++) {
                        for (int L = 0; L < bandsdim; L++) {
                            if (!(in >> aod[AOD] >> lambda[L]
                                >> Im[GG][SSA][AOD][L]
                                >> gtau[GG][SSA][AOD][L]
                                >> ag[GG][SSA][AOD][L]
                                >> btau[GG][SSA][AOD][L]
                                >> ab[GG][SSA][AOD][L])) {
                                return false;
                            }
                        }
                    }
                }
            }
            return true;
        }

        void fromFile(const std::string& filepath) {
            if (!read(filepath)) {
                std::printf("Error reading RTM aerosol LUT file: %s\n", filepath.c_str());
            }
        }

    private:
        // Allocate a 4D array [ggdim][ssadim][aoddim][bandsdim] backed by one contiguous buffer.
        MAGIC_REAL**** allocate4D(Field4D& mem) {
            const int nGG   = ggdim;
            const int nSSA  = ssadim;
            const int nAOD  = aoddim;
            const int nL    = bandsdim;

            const int nCells = nGG * nSSA * nAOD * nL;
            const int nRows  = nGG * nSSA * nAOD;    // each row is length nL
            const int nBlk3  = nGG * nSSA;           // each block3 has nAOD rows

            mem.data = new MAGIC_REAL[nCells];
            mem.p2   = new MAGIC_REAL*[nRows];
            mem.p3   = new MAGIC_REAL**[nBlk3];
            mem.p4   = new MAGIC_REAL***[nGG];

            // Build pointers
            for (int gg_i = 0; gg_i < nGG; ++gg_i) {
                mem.p4[gg_i] = mem.p3 + gg_i * nSSA;
                for (int ssa_i = 0; ssa_i < nSSA; ++ssa_i) {
                    const int blk3_index = gg_i * nSSA + ssa_i;
                    mem.p4[gg_i][ssa_i] = mem.p2 + blk3_index * nAOD;

                    for (int aod_i = 0; aod_i < nAOD; ++aod_i) {
                        const int row_index = (blk3_index * nAOD + aod_i);
                        mem.p4[gg_i][ssa_i][aod_i] = mem.data + row_index * nL;
                    }
                }
            }

            return mem.p4;
        }

        void free4D(Field4D& mem) {
            delete[] mem.p4;
            delete[] mem.p3;
            delete[] mem.p2;
            delete[] mem.data;
            mem.p4 = nullptr; mem.p3 = nullptr; mem.p2 = nullptr; mem.data = nullptr;
        }
    };

    namespace Correction {

        struct Spectral {

            // Spectral Correction of Cloud Effect
            // We must occrect from broadband cloud transmission to wavelength dependent transmission

            static constexpr int CodCols = 7;   // only uses indices 1-6

            int lmax;       // allocated rows 
            int nrows;      // number of rows read 

            // K corresponding to Cloud Optical Depth
            MAGIC_REAL k_COD[CodCols]; // [0..6], only 1..6 used

            MAGIC_REAL*  lambda;   // size lmax
            MAGIC_REAL*  corData;  // size lmax * CodCols (contiguous)
            MAGIC_REAL** cor;      // row pointers: cor[row][col]

            Spectral(int lmax_in)
                : lmax(lmax_in), nrows(0),
                lambda(nullptr), corData(nullptr), cor(nullptr) {

                // Basic alloc 
                lambda  = new MAGIC_REAL[lmax];
                corData = new MAGIC_REAL[lmax * CodCols];
                cor     = new MAGIC_REAL*[lmax];

                for (int i = 0; i < lmax; i++) {
                    cor[i] = corData + i * CodCols;
                }

                // zero this array
                for (int i = 0; i < CodCols; i++) k_COD[i] = 0.0;
            }

            ~Spectral() {
                delete[] cor;
                delete[] corData;
                delete[] lambda;
            }

            Spectral(const Spectral&) = delete;
            Spectral& operator=(const Spectral&) = delete;

            // File format:
            // line 1: header (ignored)
            // line 2: dummy + k_COD[1..6]
            // next 32 lines: lambda[kk] + cor[kk][1..6] for kk=1..32 (1-based)
            bool read(std::string& filepath) {
                std::ifstream in(filepath);
                if (!in) return false;

                // Skip header line
                std::string header;
                std::getline(in, header);

                // Read dummy + 6 k values
                MAGIC_REAL dummy;
                if (!(in >> dummy
                        >> k_COD[1] >> k_COD[2] >> k_COD[3]
                        >> k_COD[4] >> k_COD[5] >> k_COD[6])) {
                    return false;
                }

                // Read 32 rows (or as many as fit)
                nrows = 0;

                int maxRead = KATO_MAX +1;   
                if (lmax < maxRead) {
                    // Something is funky
                    // Reading more bands than there are katobands
                    return false;
                }

                for (int kk = KATO_MIN; kk < maxRead; kk++) {
                    if (!(in >> lambda[kk]
                            >> cor[kk][1] >> cor[kk][2] >> cor[kk][3]
                            >> cor[kk][4] >> cor[kk][5] >> cor[kk][6])) {
                        return false;
                    }
                    nrows++;
                }

                return true;
            }

            void fromFile(Config c) {
                
                // Get correct filepath
                std::string fp = c.path+c.lambda_cor;
                if (!read(fp)) {
                    printf("Error reading lambda correction file!\n");
                }
            }
        };

        struct Absorber {
            // Only exists for the spectral mode 
            int nrows;   // either waterdim or ozonedim
            int ncols;   // number of bands

            // 1D arrays
            MAGIC_REAL* concen;  // [nrows]
            MAGIC_REAL* lambda;  // [ncols]
            MAGIC_REAL* ag;      // [ncols]
            MAGIC_REAL* ab;      // [ncols]

            // 2D arrays (row pointers + contiguous storage)
            MAGIC_REAL** delta_global;   // [nrows][ncols]
            MAGIC_REAL** delta_beam;   // [nrows][ncols]
            MAGIC_REAL*  delta_global_data; // contiguous [nrows*ncols]
            MAGIC_REAL*  delta_beam_data; // contiguous [nrows*ncols]

            Absorber(int rows, int cols)
                : nrows(rows), ncols(cols),
                concen(nullptr), lambda(nullptr), ag(nullptr), ab(nullptr),
                delta_global(nullptr), delta_beam(nullptr),
                delta_global_data(nullptr), delta_beam_data(nullptr) {

                // 1D allocations
                concen = new MAGIC_REAL[nrows];
                lambda = new MAGIC_REAL[ncols];
                ag     = new MAGIC_REAL[ncols];
                ab     = new MAGIC_REAL[ncols];

                // 2D allocations (row pointers + contiguous block)
                delta_global = new MAGIC_REAL*[nrows];
                delta_beam = new MAGIC_REAL*[nrows];

                delta_global_data = new MAGIC_REAL[nrows * ncols];
                delta_beam_data = new MAGIC_REAL[nrows * ncols];

                for (int r = 0; r < nrows; r++) {
                    delta_global[r] = delta_global_data + r * ncols;
                    delta_beam[r] = delta_beam_data + r * ncols;
                }
            }

            ~Absorber() {
                delete[] delta_beam_data;
                delete[] delta_global_data;
                delete[] delta_beam;
                delete[] delta_global;

                delete[] ab;
                delete[] ag;
                delete[] lambda;
                delete[] concen;
            }

            // Disable copy for sanity
            Absorber(const Absorber&) = delete;
            Absorber& operator=(const Absorber&) = delete;

            // Read table:
            // - first line is ignored
            // - then nrows*ncols entries, each entry has 6 MAGIC_REALs:
            //   concen, lambda, delta_global, ag, delta_beam, ab
            bool read(const std::string& filepath) {
                std::ifstream in(filepath);
                if (!in) return false;

                // Ignore header line
                std::string header;
                std::getline(in, header);

                // Read values
                for (int r = 0; r < nrows; ++r) {
                    for (int c = 0; c < ncols; ++c) {
                        MAGIC_REAL conc, lam, dig, a_g, dib, a_b;

                        if (!(in >> conc >> lam >> dig >> a_g >> dib >> a_b)) {
                            return false; // parse error or premature EOF
                        }

                        // As in DWD code
                        concen[r]        = conc;
                        lambda[c]        = lam;
                        delta_global[r][c]   = dig;
                        ag[c]            = a_g;
                        delta_beam[r][c]   = dib;
                        ab[c]            = a_b;
                    }
                }

                return true;
            }

            void fromFile(const std::string& filepath) {
                if (!read(filepath)) {
                    std::printf("Error reading absorber table file: %s\n", filepath.c_str());
                }
            }
        };

    }

    struct LandUse {

        // Grid dimensions 
        // taken from comment in DWD code
        static constexpr unsigned int ncols = 2160;
        static constexpr unsigned int nrows = 1080;
        static constexpr unsigned int kk = ncols * nrows;

        // Raw land‑use grid
        // Each entry is one byte
        unsigned char* grid = nullptr;

        LandUse() = default;

        ~LandUse() {
            delete[] grid;
        }

        // We do not allow copies, data could be large
        LandUse(const LandUse&) = delete;
        LandUse& operator=(const LandUse&) = delete;

        // Moves OK
        LandUse(LandUse&& other) noexcept {
            *this = std::move(other);
        }

        LandUse& operator=(LandUse&& other) noexcept {
            if (this != &other) {
                grid = other.grid;
                other.grid = nullptr;
            }
            return *this;
        }

        bool read(const std::string& filepath) {

            std::ifstream in(filepath, std::ios::binary);
            if (!in) {
                return false;
            }

            // Alloc buffer if it doesn't already exist
            if (!grid) {
                grid = new unsigned char[kk];
            }

            // Read the entire file 
            in.read(reinterpret_cast<char*>(grid), kk);

            // Check expected file size
            if (in.gcount() != static_cast<std::streamsize>(kk)) {
                std::printf("ERROR: Land-use file size mismatch (%lld bytes read, expected %u)\n",
                            static_cast<long long>(in.gcount()), kk);
                return false;
            }

            return true;
        }

        void fromFile(const std::string& filepath) {
            if (!read(filepath)) {
                std::printf("Error reading land‑use climatology: %s\n",
                            filepath.c_str());
            } 
        }

        // For access via [row, col]
        inline unsigned char at(unsigned int row, unsigned int col) const {
            return grid[row * ncols + col];
        }

        // For calculating an index in the landuse file 
        inline int index(MAGIC_EXACT lat, MAGIC_EXACT lon) const {
            // Both inputs in DEGREES

            // Landuse file requires longitude [-180,180]
            MAGIC_EXACT lon_to_use;
            if (lon < 0) lon_to_use = 360. + lon;
            else lon_to_use = lon;
            
            int klu = lrint(6 * lon_to_use)
                + lrint(6 * (90 - lat + 1.0 / 12.0)) * ncols;

            return klu;

        }
        
    };



    struct Climatology {

        // Grid dims
        int nlon = 0;
        int nlat = 0;

        // Grid spacing
        MAGIC_EXACT ddeg = 0.0;

        // Axes
        MAGIC_EXACT* lon = nullptr;   // size nlon
        MAGIC_EXACT* lat = nullptr;   // size nlat

        // Field stored row-major:
        //
        //     val[ilat * nlon + ilon]
        //
        // i.e. latitude is row index, longitude is column index.
        MAGIC_REAL* val = nullptr;

        Climatology() = default;

        ~Climatology() {
            delete[] lon;
            delete[] lat;
            delete[] val;
        }

        Climatology(const Climatology&) = delete;
        Climatology& operator=(const Climatology&) = delete;

        Climatology(Climatology&& other) noexcept {
            *this = std::move(other);
        }

        Climatology& operator=(Climatology&& other) noexcept {
            if (this != &other) {
                delete[] lon;
                delete[] lat;
                delete[] val;

                lon = other.lon; other.lon = nullptr;
                lat = other.lat; other.lat = nullptr;
                val = other.val; other.val = nullptr;

                nlon = other.nlon; other.nlon = 0;
                nlat = other.nlat; other.nlat = 0;
                ddeg = other.ddeg; other.ddeg = 0.0f;
            }
            return *this;
        }

        inline MAGIC_REAL& at(int ilat, int ilon) noexcept {
            return val[ilat * nlon + ilon];
        }

        inline const MAGIC_REAL& at(int ilat, int ilon) const noexcept {
            return val[ilat * nlon + ilon];
        }

        bool read(const std::string& filepath,
                int month,
                int extra_cols_to_skip = 0) {

            if (month < 1 || month > 12) {
                std::printf("Invalid month %d, must be 1..12\n", month);
                return false;
            }

            std::ifstream in(filepath);
            if (!in) return false;

            std::string header;
            std::getline(in, header);

            MAGIC_EXACT lon0, lat0;
            MAGIC_EXACT lon_tmp, lat_tmp;
            float months[12];

            // Read first row to determine nlat.
            if (!(in >> lon0 >> lat0)) return false;

            for (int k = 0; k < 12; k++) in >> months[k];

            for (int k = 0; k < extra_cols_to_skip; k++) {
                MAGIC_EXACT dummy;
                in >> dummy;
            }

            nlat = 1;

            while (true) {
                std::streampos pos = in.tellg();

                if (!(in >> lon_tmp >> lat_tmp)) break;

                if (lon_tmp != lon0) {
                    in.seekg(pos);
                    break;
                }

                for (int k = 0; k < 12; k++) in >> months[k];

                for (int k = 0; k < extra_cols_to_skip; k++) {
                    MAGIC_EXACT dummy;
                    in >> dummy;
                }

                nlat++;
            }

            // Count total rows.
            in.clear();
            in.seekg(0);
            std::getline(in, header);

            int total_rows = 0;

            while (true) {
                if (!(in >> lon_tmp >> lat_tmp)) break;

                for (int k = 0; k < 12; k++) in >> months[k];

                for (int k = 0; k < extra_cols_to_skip; k++) {
                    MAGIC_EXACT dummy;
                    in >> dummy;
                }

                total_rows++;
            }

            nlon = total_rows / nlat;

            lon = new MAGIC_EXACT[nlon];
            lat = new MAGIC_EXACT[nlat];
            val = new MAGIC_REAL[nlon * nlat];

            // Fill actual data.
            //
            // File order is assumed to be:
            //
            //   for lon:
            //       for lat:
            //
            // We store as:
            //
            //   val[lat * nlon + lon]
            in.clear();
            in.seekg(0);
            std::getline(in, header);

            for (int ilon = 0; ilon < nlon; ilon++) {
                for (int ilat = 0; ilat < nlat; ilat++) {

                    in >> lon_tmp >> lat_tmp;

                    lon[ilon] = lon_tmp;
                    lat[ilat] = lat_tmp;

                    for (int k = 0; k < 12; k++) {
                        in >> months[k];

                        if (k == month - 1) {
                            at(ilat, ilon) = months[k];
                        }
                    }

                    for (int k = 0; k < extra_cols_to_skip; k++) {
                        MAGIC_EXACT dummy;
                        in >> dummy;
                    }
                }
            }

            // Derive the grid spacing from the axis we just read, rather than
            // assuming the longitudes tile a full 360 degrees.
            //
            // That assumption only holds for cell-centred grids. The water vapour
            // and ozone climatologies are node-centred and repeat the antimeridian
            // (-180 and +180 are both present), so nlon is one greater than the
            // number of cells and 360 / nlon understates the spacing -- by 0.07%
            // for water vapour and 0.3% for ozone, which skews every interpolation
            // weight in interpolate().
            if (nlon > 1) {
                ddeg = std::fabs(lon[nlon - 1] - lon[0])
                     / static_cast<MAGIC_EXACT>(nlon - 1);
            } else {
                ddeg = 360.0;
            }

            // interpolate() applies the same ddeg to both axes, so a file whose
            // latitude spacing differs from its longitude spacing would be
            // silently mis-weighted. None of the shipped climatologies do.
            if (nlat > 1) {
                const MAGIC_EXACT dlat = std::fabs(lat[nlat - 1] - lat[0])
                                       / static_cast<MAGIC_EXACT>(nlat - 1);

                if (std::fabs(dlat - ddeg) > 1e-6 * ddeg) {
                    std::printf("WARNING: climatology %s has unequal lat/lon spacing "
                                "(%g vs %g); interpolation assumes they match\n",
                                filepath.c_str(), dlat, ddeg);
                }
            }

            return true;
        }

        void fromFile(const std::string& filepath,
                    int month,
                    int extra_cols_to_skip = 0) {
            if (!read(filepath, month, extra_cols_to_skip)) {
                std::printf("Error reading climatology file: %s\n",
                            filepath.c_str());
            }
        }

        template <typename T>
        int borders(const T* xa, int n, T x) const {
            int klo = 0;
            int khi = n - 1;

            const T sign = (xa[n - 1] < xa[0]) ? -1.0f : 1.0f;

            while (khi - klo > 1) {
                const int k = (khi + klo) / 2;

                if (sign * xa[k] > sign * x)
                    khi = k;
                else
                    klo = k;
            }

            return khi;
        }

        MAGIC_REAL interpolate(MAGIC_EXACT lon_q, MAGIC_EXACT lat_q) const {

            const int i = borders(lon, nlon, lon_q);
            const int j = borders(lat, nlat, lat_q);

            const MAGIC_EXACT wlonl = 1. - (lon_q - lon[i - 1]) / ddeg;
            const MAGIC_EXACT wlonh = 1. - wlonl;

            const MAGIC_EXACT wlath = 1. - (lat_q - lat[j]) / ddeg;
            const MAGIC_EXACT wlatl = 1. - wlath;


            // Row-major
            //     val[ilat * nlon + ilon]
            // Overkill casting for safety
            return static_cast<MAGIC_REAL>( wlath * wlonl  * at(j,     i - 1) +
                ( wlatl * wlonl ) * at(j - 1, i - 1) +
                ( wlath * wlonh ) * at(j,     i    ) +
                ( wlatl * wlonh ) * at(j - 1, i    ) );
        }
    };


    struct AerosolOptics {

        // The existing aerosol climatology has AOD varying with lat/lon
        // but NOT these. These are (assumed) constant in space 
        // We do not want to store a huge table of constants
        // So we do this instead

        MAGIC_REAL ssa = 0.95;   // single scattering albedo
        MAGIC_REAL gg = 0.7;     // asymmetry param

    };

    struct Aerosol {

        Climatology aod; 
        AerosolOptics optics;
    };

    // Only exists to HOLD the existing climatologies 
    // So we don't mess up variable names
    // there are many variables named water, ozone, etc...
    struct Climate {

        Climatology water;
        Climatology ozone;
        Aerosol aero;
        LandUse lu;

        // Reads in from file...
        Climate(Config c, int month) {

            // Ozone climatology
            ozone.fromFile(c.path+c.o3clim, month);

            // Water vapour climatology
            water.fromFile(c.path+c.hclim, month);

            // Aerosol climatology
            aero.aod.read(c.path+c.aeroclim, month, 12);

            // Landuse climatology
            lu.fromFile(c.path+c.landuse_image);


        }

        void makeLocalClimate(PixelClimate& clim, Area a, bool talk = false) {


            clim.aod = this->aero.aod.interpolate(a.degrees_lon, a.degrees_lat);
            if (talk) printf("Got aod %g \n", clim.aod);

            clim.ssa = this->aero.optics.ssa;
            clim.gg = this->aero.optics.gg;
            if (talk) printf("Got ssa %g \n", clim.ssa);

            clim.h2o = this->water.interpolate(a.degrees_lon, a.degrees_lat);
            if (talk) printf("Got water %g \n", clim.h2o);

            // PASSING DEGREES_LON AND DEGREES_LAT SHOULD BE CORRECT PHYSICS
            // TODO: Check
            clim.o3 = this->ozone.interpolate(a.degrees_lon, a.degrees_lat);
            if (talk) printf("Got ozone %g with lon %f and lat %f \n", clim.o3, a.degrees_lon, a.degrees_lat);

        }
            


    };



}