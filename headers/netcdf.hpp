#pragma once 
#include <cstdlib>
#include <vector>

#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

#define NC_CHECK(call) do {                           \
	int status = (call);                                \
	if (status != NC_NOERR) {                            \
		fprintf(stderr, "%s:%d: %s\n",                     \
			__FILE__, __LINE__, nc_strerror(status));  \
			return status;                                     \
		}                                                    \
	} while(0)

template <typename T>
int nc_put_vara_magic_int(int ncid, int varid, size_t* start,
    size_t* count, T* data) {

    if constexpr (std::is_same_v<T, short>) {
        return nc_put_vara_short(ncid, varid, start, count, data);
    }
    else if constexpr (std::is_same_v<T, int>) {
        return nc_put_vara_int(ncid, varid, start, count, data);
    }
    else {
        static_assert(std::is_integral_v<T>,
                      "Unsupported MAGIC_INT type for NetCDF output");
    }
}

template <typename T>
int nc_put_att_magic_int(int ncid, int varid, const char* name,
                         const T* value) {

    if constexpr (std::is_same_v<T, short>) {
        return nc_put_att_short(ncid, varid, name, NC_SHORT, 1, value);
    }
    else if constexpr (std::is_same_v<T, int>) {
        return nc_put_att_int(ncid, varid, name, NC_INT, 1, value);
    }
    else {
        static_assert(std::is_integral_v<T>,
                      "Unsupported MAGIC_INT type for NetCDF attribute");
    }
}
template <typename T>
int nc_put_var_magic_real(int ncid, int varid, const T* data, std::size_t n) {
    
    if constexpr (std::is_same_v<T, float>) {
        return nc_put_var_float(ncid, varid, data);
    }
    else if constexpr (std::is_same_v<T, double>) {
        return nc_put_var_double(ncid, varid, data);
    }
    else if constexpr (std::is_same_v<T, long double>) {
        std::vector<double> tmp(n);

        for (std::size_t i = 0; i < n; ++i) {
            tmp[i] = static_cast<double>(data[i]);
        }

        return nc_put_var_double(ncid, varid, tmp.data());
    }
    else {
        static_assert(std::is_floating_point_v<T>,
                      "Unsupported floating-point type for NetCDF output");
    }
}

	