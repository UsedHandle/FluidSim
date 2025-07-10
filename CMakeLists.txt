cmake_minimum_required(VERSION 3.24)

set(EXE_NAME "FluidSim")

project(${EXE_NAME} LANGUAGES C CXX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

file(
	DOWNLOAD
	https://github.com/cpm-cmake/CPM.cmake/releases/download/v0.40.8/CPM.cmake
	${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake
	EXPECTED_HASH SHA256=78ba32abdf798bc616bab7c73aac32a17bbd7b06ad9e26a6add69de8f3ae4791
)
include(${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake)

set(LIBS)
CPMAddPackage("gh:libsdl-org/SDL#8d60435")
list(APPEND LIBS SDL3::SDL3)

add_executable(${EXE_NAME} src/main.cpp src/mesh.cpp src/fvsolver.cpp src/triangle.cpp vendor/imgui/imgui.cpp vendor/imgui/imgui_demo.cpp vendor/imgui/imgui_draw.cpp vendor/imgui/imgui_tables.cpp vendor/imgui/imgui_widgets.cpp vendor/imgui/imgui_impl_sdl3.cpp vendor/imgui/imgui_impl_sdlrenderer3.cpp)
target_include_directories(${EXE_NAME} PRIVATE vendor/ vendor/imgui/ )
target_compile_features(${EXE_NAME} PRIVATE cxx_std_23)
target_link_libraries(${EXE_NAME} PRIVATE ${LIBS})
