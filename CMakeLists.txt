cmake_minimum_required(VERSION 3.18)
project(assignment_tests VERSION 0.1 LANGUAGES C)

# write compile_commands.json file
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

# Set a default build type if none was specified
set(default_build_type "Debug")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set(SRC
    src/main.c
    )

add_executable(assignment_tests ${SRC})
set_property(TARGET assignment_tests PROPERTY C_STANDARD 11)

# compile flags
target_compile_options(
    assignment_tests PRIVATE 
    $<$<CONFIG:RELEASE>:-O3 -march=native>
    $<$<CONFIG:DEBUG>:-O0 -Wall -Werror>
    )

# 3rd party
find_package(HDF5 REQUIRED COMPONENTS C HL)
target_include_directories(assignment_tests PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(assignment_tests PRIVATE ${HDF5_C_LIBRARIES} ${HDF5_C_HL_LIBRARIES})
target_compile_definitions(assignment_tests PRIVATE ${HDF5_DEFINITIONS})

find_package(GSL REQUIRED)
target_link_libraries(assignment_tests PRIVATE GSL::gsl GSL::gslcblas)

