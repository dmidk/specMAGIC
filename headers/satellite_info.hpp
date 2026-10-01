#pragma once 
#include "toml.hpp"
#include "navigation.hpp"
#include <filesystem>
#include <cstdlib>
#include <stdexcept>
#include <unordered_map>


struct Metadata {

    // number of pixels in x and y
    unsigned int num_lines;
    unsigned int num_columns;

    // Defines the origin of the image
    int line_offset;
    int column_offset;

    // get this from the image
    int resolution; // in metres!

    MAGIC_EXACT satellite_radius_km;
    MAGIC_EXACT subsatellite_longitude_deg;
    MAGIC_EXACT full_disk_scan_min;

    // Selected sampling angle for this run
    MAGIC_EXACT angular_sampling_rad;

    int wavelength;

    bool flip_vertical;
    bool flip_horizontal;

    // Available sampling grids for this satellite
    std::unordered_map<int, MAGIC_EXACT> angular_sampling;

    // available wavelengths for this satellite
    std::unordered_map<std::string, int> wavelengths;
};

inline Metadata loadMetadata(const std::string& channel){
    const char* satellite_env = std::getenv("SATELLITE");

    if (satellite_env == nullptr || std::string(satellite_env).empty()) {
        throw std::runtime_error(
            "SATELLITE environment variable is not set");
    }

    // look for the toml
    std::string satellite = satellite_env;
    std::filesystem::path config_path = std::filesystem::path("satellites") / (satellite + ".toml");

    // there has to be a toml
    if (!std::filesystem::exists(config_path)) {
        throw std::runtime_error(
            "Satellite configuration file not found: " +
            config_path.string());
    }

    // read the toml
    toml::table config = toml::parse_file(config_path.string());

    Metadata metadata;

    // put the toml values into the metadata struct
    metadata.satellite_radius_km = config["navigation"]["satellite_radius_km"].value_or(0.0);
    {
        auto parsed_lon = config["navigation"]["subsatellite_longitude_deg"].value<MAGIC_EXACT>();
        if (!parsed_lon) {
            throw std::runtime_error(
                "Missing required key 'navigation.subsatellite_longitude_deg' in " +
                config_path.string());
        }
        metadata.subsatellite_longitude_deg = *parsed_lon;
    }
    metadata.full_disk_scan_min = config["navigation"]["full_disk_scan_min"].value_or(0.0);
    
    metadata.flip_horizontal = config["orientation"]["flip_horizontal"].value_or(false);
    metadata.flip_vertical = config["orientation"]["flip_vertical"].value_or(false);

    // to deal with the various possible resolutions
    if (auto* grids = config["grids"].as_table()) {
        for (auto& [key, value] : *grids) {
            int resolution = std::stoi(std::string(key));

            if (auto* grid = value.as_table()) {
                if (auto sampling = (*grid)["angular_sampling_rad"].value<MAGIC_EXACT>()) {
                    metadata.angular_sampling[resolution] = *sampling;
                }
            }
        }
    }

    // to deal with the various possible channels
    if (auto* channels = config["channels"].as_table()) {
        for (auto& [key, value] : *channels) {
            if (auto wavelength = value.value<int>()) {metadata.wavelengths[
                    std::string(key)] = *wavelength;
            }
        }
    }

    // what channel are we using now
    auto wavelength_it = metadata.wavelengths.find(channel);

    if (wavelength_it == metadata.wavelengths.end()) {
        throw std::runtime_error("Channel '" + channel + "' is not defined in " + config_path.string());
    }

    metadata.wavelength = wavelength_it->second;

    return metadata;
}