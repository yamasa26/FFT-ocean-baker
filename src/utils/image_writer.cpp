#include "image_writer.hpp"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <vector>
#include <algorithm>
#include <filesystem>

namespace fs = std::filesystem;

void save_heightmap(const std::string& filename, int width, int height, const std::vector<float>& data) {
    // 保存先の親ディレクトリが存在するか確認して、なければ作成
    fs::path filepath(filename);
    if (filepath.has_parent_path()) {
        fs::create_directories(filepath.parent_path());
    }

    std::vector<unsigned char> pixels(width * height);
    
    // データの最小値と最大値を見つける
    auto [min_it, max_it] = std::minmax_element(data.begin(), data.end());
    float min_val = *min_it;
    float max_val = *max_it;
    float range = max_val - min_val;

    for (int i = 0; i < width * height; ++i) {
        // 0.0 ~ 1.0 に正規化してから 0 ~ 255 に変換
        float normalized = (range > 0.0001f) ? (data[i] - min_val) / range : 0.0f;
        pixels[i] = static_cast<unsigned char>(normalized * 255.0f);
    }

    stbi_write_png(filename.c_str(), width, height, 1, pixels.data(), width);
}

void save_image_rgb(const std::string& filename, int width, int height, 
                    const std::vector<float>& r, 
                    const std::vector<float>& g, 
                    const std::vector<float>& b) {
    fs::path filepath(filename);
    if (filepath.has_parent_path()) fs::create_directories(filepath.parent_path());

    std::vector<unsigned char> pixels(width * height * 3);
    for (int i = 0; i < width * height; ++i) {
        // -1.0 ~ 1.0 の範囲を 0 ~ 255 に変換（法線用）
        pixels[i * 3 + 0] = static_cast<unsigned char>((r[i] * 0.5f + 0.5f) * 255.0f);
        pixels[i * 3 + 1] = static_cast<unsigned char>((g[i] * 0.5f + 0.5f) * 255.0f);
        pixels[i * 3 + 2] = static_cast<unsigned char>((b[i] * 0.5f + 0.5f) * 255.0f);
    }
    stbi_write_png(filename.c_str(), width, height, 3, pixels.data(), width * 3);
}

void save_combined_rgba(const std::string& filename, int width, int height, 
                        const std::vector<float>& h, 
                        const std::vector<float>& nx, 
                        const std::vector<float>& ny, 
                        const std::vector<float>& foam) {
    fs::path filepath(filename);
    if (filepath.has_parent_path()) fs::create_directories(filepath.parent_path());

    std::vector<unsigned char> pixels(width * height * 4);
    
    // Heightの正規化用
    const float h_min = -5.0f; 
    const float h_max = 5.0f;

    for (int i = 0; i < width * height; ++i) {
        // R: Height (0~1に正規化)
        pixels[i * 4 + 0] = static_cast<unsigned char>(std::clamp((h[i] - h_min) / (h_max - h_min), 0.0f, 1.0f) * 255.0f);
        // G: Normal X (-1~1 -> 0~1)
        pixels[i * 4 + 1] = static_cast<unsigned char>((nx[i] * 0.5f + 0.5f) * 255.0f);
        // B: Normal Y (-1~1 -> 0~1)
        pixels[i * 4 + 2] = static_cast<unsigned char>((ny[i] * 0.5f + 0.5f) * 255.0f);
        // A: Foam (0~1)
        pixels[i * 4 + 3] = static_cast<unsigned char>(std::clamp(foam[i], 0.0f, 1.0f) * 255.0f);
    }
    stbi_write_png(filename.c_str(), width, height, 4, pixels.data(), width * 4);
}
