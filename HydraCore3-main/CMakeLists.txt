cmake_minimum_required (VERSION 3.10)
project (test)

#set(CMAKE_BUILD_TYPE Release)
# get rid of annoying MSVC warnings.
add_definitions(-D_CRT_SECURE_NO_WARNINGS)

set (CMAKE_CXX_STANDARD 17)

include_directories(".")
if(WIN32)
  link_directories(${ADDITIONAL_LIBRARY_DIRS})
endif()

#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}") # to find 'FindFreeImage.cmake' script
#find_package(FreeImage REQUIRED)

find_package(OpenMP)
find_package(Vulkan)

#include_directories(${Vulkan_INCLUDE_DIR})
#include_directories(${FREEIMAGE_INCLUDE_DIRS})
include_directories(#"../vkutils"
                    #"../volk" 
                    "../LiteMath" 
                    "../loader_utils"
                    "../geo_utils"
                    "../RayTracing"
                    "../utils")
                    
#link_directories("../volk")
link_directories("../RayTracing/lib")

#add_compile_definitions(USE_VOLK)

#uncomment this to detect broken memory problems via gcc sanitizers
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fsanitize-address-use-after-scope -fno-omit-frame-pointer -fsanitize=leak -fsanitize=undefined -fsanitize=bounds-strict")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -mavx2 -fPIC")

add_executable(hydra main.cpp 
                     integrator_pt1.cpp
                     integrator_pt2.cpp
                     integrator_pt_scene.cpp
                     #imageutils.cpp 
                     Bitmap.cpp
                     ../LiteMath/texture2d.cpp
                     #integrator_pt1_generated.cpp
                     #integrator_pt1_generated_init.cpp 
                     #integrator_pt1_generated_ds.cpp 
                     ../loader_utils/pugixml.cpp
                     ../loader_utils/hydraxml.cpp
                     #../loader_utils/scene_mgr.cpp
                     ../geo_utils/cmesh.cpp
                     addsrc/DSTree/dstree.cpp
                     #../RayTracing/EmbreeRT.cpp #comment either this or dstree.cpp
                     )

set_target_properties(hydra PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
target_link_libraries(hydra LINK_PUBLIC  OpenMP::OpenMP_CXX  #dl     ${FREEIMAGE_LIBRARIES}
                      #embree3 embree_sse42 embree_avx embree_avx2 lexers simd sys tbb # embree
                      ) 
