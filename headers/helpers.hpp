#pragma once 
#include <string>


struct Config {

    std::string aeroclim, hclim, o3clim;
    std::string landuse_image, clut_spec, clut_spec_h2o, clut_spec_o3;
    std::string lambda_cor, cloud_list;
    std::string out_global_path, out_beam_path, out_cal_path, out_clear_path;
    std::string gr_alb_file;
    std::string modis_brdf_dir;
    std::string path;

    int xadim{}, yadim{}, xhdim{}, yhdim{}, xo3dim{}, yo3dim{};
    int ggdim{}, ssadim{}, aoddim{}, h2odim{}, o3dim{}, num_bands{};
    int latdim{}, londim{};
    float latbegin{}, lonbegin{}, dxy{}, deltalon{}, iconres{};
    int iconflag{};
    int use_modis_brdf_albedo{0};

 
};

inline void printPixels(size_t day, size_t night, size_t invalid) {

    bool concern = false;

    float perc1 = static_cast<float>(invalid) / (day + night + invalid);
    float perc2 = static_cast<float>(day) / (day + night + invalid);
    float perc3 = static_cast<float>(night) / (day + night + invalid);

    if (perc1 * 100 > 6) concern = true;

    if (concern) printf("[WARNING]: Many invalid pixels. Results may not be sound. \n");

    printf("Found %.1f%% of pixels to be daytime, %.1f%% to be nighttime, and %.1f%% to be invalid. \n", perc2 * 100, 
        perc3 * 100, perc1 * 100);
    
}

inline int channelToWavelength(std::string& channel) {
    
    if (channel == "vis_04") return 440;
    else if (channel == "vis_05") return 510;
    else if (channel == "vis_06") return 640;
    else if (channel == "vis_08") return 865;
    else if (channel == "vis_09") return 914;
    else if (channel == "nir_13") return 1380;
    else if (channel == "nir_16") return 1610;
    else if (channel == "nir_22") return 2250;
    else {
        printf("Invalid satellite channel!");
        return -1; 
    }
}