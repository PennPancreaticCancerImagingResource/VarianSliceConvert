cmake_minimum_required(VERSION 3.10.2)

project(VarianSliceConvert)

# Find ITK.
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable(VarianSliceConvert VarianSliceConvert.cxx )

target_link_libraries(VarianSliceConvert ${ITK_LIBRARIES})
