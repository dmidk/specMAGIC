#include <stdio.h>
#include <omp.h>
#include <string>
#include <cstdlib>
#include <stdexcept>
#include <fstream>
#include "../headers/types.hpp"
#include "../headers/read.hpp"
#include "../headers/image.hpp"
#include "../headers/helpers.hpp"
#include "../headers/run.hpp"


int main(int argc, char* argv[]) {

    // Start the timer
    double start_time = omp_get_wtime();

    int threads;

    // get the thread count from environment
    #pragma omp parallel
    {
        threads = omp_get_num_threads();
    }

    // Safety check
    if (argc != 4) {
        fprintf(stderr,
            "Usage: %s <home-directory> <channel> <timer-flag>\n"
            "Expected exactly 3 input arguments, but received %d.\n"
            "Example: %s /path/to/home ir_105 1\n",
            argv[0], argc - 1, argv[0]);
        return 1;
    }

    // Home directory same as in the driver script
    std::string home = argv[1];
    std::string channel = argv[2];     // the mtg channel
    bool timer = std::stoi(argv[3]) != 0;

    // Get setup from the configuration file
    Config c = loadConfig("magic-config.asc", home);

    // Read in the actual satellite image and its metadata
    Image img;
    img.initfromFile(c.path+c.cloud_list);
    img.readImage(channel, satelliteFromEnv(), home);

    AlbedoType albedo_source = parseAlbedoSource();

    // Optional input fields on the output grid, one NetCDF file per table
    InputFields::Paths paths;
    paths.aerosol = envString("AEROSOL_FILE");
    paths.ozone = envString("OZONE_FILE");
    paths.water_vapour = envString("WATER_VAPOUR_FILE");
    paths.albedo = envString("ALBEDO_FILE");

    if (!paths.albedo.empty()) {
        if (albedo_source == AlbedoType::MODIS) {
            throw std::runtime_error("ALBEDO=MODIS and ALBEDO_FILE cannot be used together");
        }
        albedo_source = AlbedoType::INPUT;
    }

    Geography geo = outputGrid(c, img.info.subsatellite_longitude_deg);

    InputFields::Inputs inputs;
    inputs.request.strict = envFlag("STRICT_INPUTS");
    inputs.load(paths, geo);
    inputs.check();

    MagicRun run = runMagic(c, img, geo, albedo_source, inputs);

    // Output information to file, recording which input fields were used
    run.radiation.toFile(img.timestamp, geo, c, run.attributes);

    // Stop the timer
    double end_time = omp_get_wtime();

    // Only writeout a time if we ask for it
    if (timer) {

        std::ofstream tdat(home+"time.txt", std::ios_base::app);
        tdat << threads << " " << run.loop_seconds << " " << end_time - start_time << std::endl;
        tdat.close();

    }

    printf("Magic finished! Took %.1f seconds. \n", end_time - start_time);



    return 0;

}
