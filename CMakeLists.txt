cmake_minimum_required(VERSION 3.10)
project(FFTOceanBaker)

set(CMAKE_CXX_STANDARD 17)
add_definitions(-D_CRT_SECURE_NO_WARNINGS)

#相対パス指定で原因不明のビルドエラーが出るため、絶対パス指定で実装
#絶対パス指定は非推奨。原因がわかり次第変更する
set(VCPKG_INCLUDE_DIR "C:/vcpkg/installed/x64-windows/include")
set(VCPKG_LIB_DIR     "C:/vcpkg/installed/x64-windows/lib")

add_executable(ocean_baker 
    src/main.cpp 
    src/ocean/simulator.cpp
    src/utils/image_writer.cpp
)

target_include_directories(ocean_baker PRIVATE 
    "${VCPKG_INCLUDE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/ocean"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/utils"
)

target_link_libraries(ocean_baker PRIVATE 
    "${VCPKG_LIB_DIR}/fftw3.lib"
    "${VCPKG_LIB_DIR}/fftw3f.lib"
)