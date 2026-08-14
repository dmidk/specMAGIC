#include "../headers/write.hpp"

std::string makeFilename(DateTime& timestamp,  std::string& variable_name, std::string& out_path) {
 
    std::ostringstream oss;
    oss << out_path << "/"
        << variable_name << timestamp.year
        << std::setw(2) << std::setfill('0') << timestamp.month
        << std::setw(2) << std::setfill('0') << timestamp.day
        << std::setw(2) << std::setfill('0') << timestamp.hour
        << std::setw(2) << std::setfill('0') << timestamp.minute
        << ".nc";

    std::string fname = oss.str();

    return fname;

   

}

int writeToNetcdf(const std::string& filename,
    Geography& geo, Matrix& box, int year, int month, int day,
    double seconds_since_midnight, std::string name) {
    int ncid;

    // Create file
    NC_CHECK(nc_create(filename.c_str(), NC_NETCDF4 | NC_CLOBBER, &ncid));

    // Dimensions
    int time_dim, lat_dim, lon_dim;
    NC_CHECK(nc_def_dim(ncid, "time", NC_UNLIMITED, &time_dim));
    NC_CHECK(nc_def_dim(ncid, "lat",  geo.nlat, &lat_dim));
    NC_CHECK(nc_def_dim(ncid, "lon",  geo.nlon, &lon_dim));

    // Variables
    int time_var, lat_var, lon_var, matrix_var;
    int dims_time[1]      = {time_dim};
    int dims_lat[1]       = {lat_dim};
    int dims_lon[1]       = {lon_dim};
    int dims_matrix[3]    = {time_dim, lat_dim, lon_dim};

    NC_CHECK(nc_def_var(ncid, "time", NC_DOUBLE, 1, dims_time, &time_var));
    NC_CHECK(nc_def_var(ncid, "lat", NC_DOUBLE, 1, dims_lat, &lat_var));
    NC_CHECK(nc_def_var(ncid, "lon", NC_DOUBLE, 1, dims_lon, &lon_var));


    NC_CHECK(nc_def_var(ncid, name.c_str(), NC_SHORT, 3, dims_matrix, &matrix_var));

    // Attributes
    char time_units[64];
    std::snprintf(time_units, sizeof(time_units),
        "seconds since %04d-%02d-%02d 00:00:00",
        year, month, day);

    NC_CHECK(nc_put_att_text(ncid, time_var, "units",
                             std::strlen(time_units), time_units));
    NC_CHECK(nc_put_att_text(ncid, lat_var, "units", 13, "degrees_north"));
    NC_CHECK(nc_put_att_text(ncid, lon_var, "units", 12, "defgrees_east"));
    NC_CHECK(nc_put_att_text(ncid, matrix_var, "units", 6, "W m-2"));

    MAGIC_INT fill = static_cast<MAGIC_INT>(-1);
    NC_CHECK(nc_put_att_magic_int(ncid, matrix_var, "_FillValue", &fill));

    // End define mode
    NC_CHECK(nc_enddef(ncid));

     // Write coordinates
    NC_CHECK(nc_put_var_magic_real(ncid, lat_var, geo.lat,
                                static_cast<std::size_t>(geo.nlat)));

    NC_CHECK(nc_put_var_magic_real(ncid, lon_var, geo.lon,
                                static_cast<std::size_t>(geo.nlon)));

    // Write time
    size_t t0 = 0;
    NC_CHECK(nc_put_var1_double(ncid, time_var, &t0, &seconds_since_midnight));

    // Write slab
    size_t start[3] = {0, 0, 0};
    size_t count[3] = {1,
                       static_cast<size_t>(geo.nlat),
                       static_cast<size_t>(geo.nlon)};

    NC_CHECK(nc_put_vara_magic_int(ncid, matrix_var, start, count,
        box.data()));


    nc_close(ncid);

    return 0;
} 