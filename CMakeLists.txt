cmake_minimum_required(VERSION 3.11)
project(Sharp_features_3D_image)

set(CMAKE_CXX_STANDARD 11)

add_definitions(-DCGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
                -DCGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS)

if(CGAL_MESH_3_VERBOSE)
  add_definitions(-DCGAL_MESH_3_VERBOSE)
endif()

find_package(CGAL REQUIRED COMPONENTS ImageIO)

option(CGAL_ACTIVATE_CONCURRENT_MESH_3 "Activate parallelism in Mesh_3" OFF)

if(CGAL_ACTIVATE_CONCURRENT_MESH_3 OR "$ENV{CGAL_ACTIVATE_CONCURRENT_MESH_3}")
  add_definitions(-DCGAL_CONCURRENT_MESH_3)
  find_package(TBB REQUIRED)
  include(CGAL_TBB_support)
endif()

# Use Eigen
find_package(Eigen3 3.1.0 QUIET) #(3.1.0 or greater)
include(CGAL_Eigen3_support)
if(NOT TARGET CGAL::Eigen3_support)
  message("NOTICE: All examples require the Eigen3 library, and will not be compiled.")
  return()
endif()

# Use DGtal
find_package(DGtal REQUIRED)
include_directories(${DGTAL_INCLUDE_DIRS})
link_directories(${DGTAL_LIBRARY_DIRS})

# Use Boost
find_package(Boost REQUIRED COMPONENTS system filesystem program_options)

set(LINK_LIBRARIES ${DGTAL_LIBRARIES} CGAL::Eigen3_support Boost::filesystem Boost::program_options)

if(TARGET CGAL::CGAL_ImageIO)
  if(CGAL_ImageIO_USE_ZLIB)

    create_single_source_cgal_program("extract_sharp_features.cpp")
    target_link_libraries(extract_sharp_features PUBLIC ${LINK_LIBRARIES} Boost::program_options)
    target_sources(extract_sharp_features
      PRIVATE
      helper/dgtal_helper_types.h
      helper/dgtal_helper.h
      helper/timing_helper.h
      helper/image_helper.h
      helper/cgal_helper.h
      helper/geometry_helper.h
      helper/mesh_generator.h
      surface/Surfel.h
      surface/Surface_process.h
      feature_graph/surfel_selector.h
      feature_graph/polyline_extractor.h
      feature_graph/polyline_regularizer.h
      feature_graph/polyline_optimizer.h
      feature_graph/Surfel_point.h
      feature_graph/energy/Get_plane_intersection_energy.h
      feature_graph/energy/Get_smoothing_and_plane_intersection_energy.h
      element_graph/graph_element.h
      element_graph/Element_graph.h
      element_graph/Element_graph_process.h
      sharpness_measure.h
      polyline/Polylines.h
      polyline/polyline_process.h
      polyline/polyline_snap.h
      accuracy_measure/accuracy_measure.h
      sharp_feature_extraction.h
    )

    if(MSVC)
      target_compile_options(extract_sharp_features PRIVATE /W3)
    else()
      target_compile_options(extract_sharp_features PRIVATE -Wall -Wextra -Wpedantic -Werror)
    endif()

  else()
    message(STATUS "NOTICE: The example extract_sharp_features.cpp need CGAL_ImageIO to be configured with ZLIB support, and will not be compiled.")
  endif()
else()
  message(STATUS "NOTICE: Some examples need the CGAL_ImageIO library, and will not be compiled.")
endif()

if(CGAL_ACTIVATE_CONCURRENT_MESH_3 AND TARGET CGAL::TBB_support)
  foreach(
    target
    mesh_construct_initial_points)
    if(TARGET ${target})
      target_link_libraries(${target} PUBLIC CGAL::TBB_support)
    endif()
  endforeach()
endif()
