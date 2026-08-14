#pragma once
#include <string>
#include <fstream>
#include <cassert>
#include <cctype>
#include <cstring>   // std::strncpy, std::strlen
#include <cstdlib>   // std::atoi
#include <netcdf.h>
#include "matrix.hpp"
#include "netcdf.hpp"
#include "helpers.hpp"

struct Metadata {
    
    unsigned int num_lines;       // number of lines
    unsigned int num_columns;     // number of columns with information
    int line_offset;
    int column_offset;		       
    unsigned int nav_lres;	       // line resolution
    unsigned int nav_cres;	       // column resolution

    int wavelength;       // what wavelength satellite sees in

    // Default ctor so compiler doesn't complain
    Metadata() : num_lines(11136), num_columns(11136), line_offset(5568), column_offset(5568), nav_lres(222), nav_cres(222), wavelength(640) {}

    // These values are taken from DWD
    Metadata(float lambda) : num_lines(11136), num_columns(11136), line_offset(5568), column_offset(5568), nav_lres(222), nav_cres(222), wavelength(lambda) {}


};

struct Image {

    Metadata info;

    // Path to the image
    std::string image_name;

    // year, month, day, hour, minute
    DateTime timestamp;

    // False if not explicitly given 
    bool has_timestamp = false;

    // The actual image itself
    Matrix im;

    // Should be called as part of constructor
    void set_timestamp(const DateTime& dt) {
        timestamp = dt;
        has_timestamp = true;
    }

    // Reads one filename from a list file and parses the start timestamp
    bool read_from_list(const std::string& list_path) {
        std::ifstream in(list_path);
        if (!in)
            return false;

        // Read first non-empty line
        while (std::getline(in, image_name)) {
            if (!image_name.empty())
                break;
        }

        if (image_name.empty())
            return false;

        // Expected something like:
        // fill-vis_06-20260504114000-20260504115000-71.nc
        //
        // Parse start timestamp: YYYYMMDDHHMMSS
        std::size_t last_slash = image_name.find_last_of("/\\");
        const std::string& name =
            (last_slash == std::string::npos)
                ? image_name
                : image_name.substr(last_slash + 1);

        std::size_t first_dash  = name.find('-');
        std::size_t second_dash = name.find('-', first_dash + 1);

        if (second_dash == std::string::npos)
            return false;

        std::size_t ts_start = second_dash + 1;

        if (ts_start + 14 > name.size())
            return false;

        DateTime dt;
        if (!parse_datetime(name.substr(ts_start, 14), dt))
            return false;

        set_timestamp(dt);
        return true;
    }

    void initfromFile(const std::string& list_path) {

        // Initialises date time info from a file path
        // the actual image is alloced and filled in later

        bool success = read_from_list(list_path);
        if (!success) printf("Could not read image location information! \n");

        return;

    }

    void readImage(std::string channel) {

        int lambda = channelToWavelength(channel);
        // Initialise the metadata
        info = Metadata(lambda);

        // Read the actual image
        int status = parseImage(channel); 

        assert(status < 1 && "Failed to read in satellite image!");

    }

private:
    // Parses YYYYMMDDHHMM (seconds ignored if present)
    static inline bool parse_datetime(
        const std::string& s,
        DateTime& dt
    ) {
        if (s.size() < 12)
            return false;

        auto to_uint = [](const std::string& v) -> unsigned int {
            unsigned int result = 0;
            for (char c : v) {
                if (!std::isdigit(static_cast<unsigned char>(c)))
                    return 0;
                result = result * 10 + (c - '0');
            }
            return result;
        };

        dt.year   = to_uint(s.substr(0, 4));
        dt.month  = to_uint(s.substr(4, 2));
        dt.day    = to_uint(s.substr(6, 2));
        dt.hour   = to_uint(s.substr(8, 2));
        dt.minute = to_uint(s.substr(10, 2));

        return dt.valid();
    }

    int parseImage(std::string channel) {

        // Note: in general, we avoid size_t for math purposes, because its bytesize is platform dependent. 
        // However the netcdf library expects it. So we use it here.

        // For getting the netcdf info
        int ncid = -1, dimid = -1, varid = -1;
        
        NC_CHECK(nc_open(image_name.c_str(), NC_NOWRITE, &ncid));

        // Read x, y dimension length
        size_t nx = -1, ny = -1;
        NC_CHECK(nc_inq_dimid(ncid, "x", &dimid));
        NC_CHECK(nc_inq_dimlen(ncid, dimid, &nx));

        NC_CHECK(nc_inq_dimid(ncid, "y", &dimid));
        NC_CHECK(nc_inq_dimlen(ncid, dimid, &ny));

        im = Matrix((int)ny, (int)nx);

        // Get the variable id
        NC_CHECK(nc_inq_varid(ncid, channel.c_str(), &varid));

        // Read sat data into short buffer
        NC_CHECK(nc_get_var_short(ncid, varid, im.data()));

        // Buffer to hold time
        char start_time[64] = {0};

        // First, inquire attribute type (string vs char)
        nc_type atype;
        size_t alen;
        NC_CHECK(nc_inq_atttype(ncid, varid, "start_time", &atype));
        NC_CHECK(nc_inq_attlen(ncid, varid, "start_time", &alen));

        if (atype == NC_STRING) {
            // NetCDF allocates memory for NC_STRING attributes
            char* tmp = nullptr;
            NC_CHECK(nc_get_att_string(ncid, varid, "start_time", &tmp));

            if (tmp) {
                std::strncpy(start_time, tmp, sizeof(start_time) - 1);
                start_time[sizeof(start_time) - 1] = '\0';
            }

            NC_CHECK(nc_free_string(1, &tmp));
        }
        else if (atype == NC_CHAR) {
            // Fixed-length character attribute
            const size_t n = std::min(alen, sizeof(start_time) - 1);
            NC_CHECK(nc_get_att_text(ncid, varid, "start_time", start_time));
            start_time[n] = '\0';
        }
        else {
            // Unexpected attribute type
            start_time[0] = '\0';
        }
        if (start_time[0] != '\0') {
            // Expected format: YYYY-MM-DD[T ]HH:MM...

            if (std::strlen(start_time) >= 16 &&
                start_time[4] == '-' &&
                start_time[7] == '-' &&
                (start_time[10] == 'T' || start_time[10] == ' ') &&
                start_time[13] == ':') {
                const int year   = std::atoi(start_time + 0);
                const int month  = std::atoi(start_time + 5);
                const int day    = std::atoi(start_time + 8);
                const int hour   = std::atoi(start_time + 11);
                const int minute = std::atoi(start_time + 14);

                timestamp = DateTime(year, month, day, hour, minute);
                has_timestamp = true;
            }
            else {
                printf("Warning: unrecognised start_time format: %s\n", start_time);
                has_timestamp = false;
            }
        } else {
            printf("Warning: start_time attribute missing or empty\n");
            has_timestamp = false;
        }

        nc_close(ncid);


        return 0;

    }
};
