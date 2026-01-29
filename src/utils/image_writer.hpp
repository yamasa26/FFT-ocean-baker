#pragma once
#include <string>
#include <vector>

void save_heightmap(
    const std::string& filename, 
    int width, 
    int height, 
    const std::vector<float>& data
);

void save_image_rgb(const std::string& filename, int width, int height, 
                    const std::vector<float>& r, 
                    const std::vector<float>& g, 
                    const std::vector<float>& b
);

void save_combined_rgba(const std::string& filename, int width, int height, 
                        const std::vector<float>& h, 
                        const std::vector<float>& nx, 
                        const std::vector<float>& ny, 
                        const std::vector<float>& foam);