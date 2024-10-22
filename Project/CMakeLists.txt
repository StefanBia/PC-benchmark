cmake_minimum_required(VERSION 3.29)
project(Project)

set(CMAKE_CXX_STANDARD 20)

# Include directories for ImGui and ImPlot
include_directories(
        imgui-master
        imgui-master/backends
        implot-master/implot
)

# Add source files for ImGui and ImPlot
file(GLOB IMGUI_SOURCES
        imgui-master/*.cpp
        imgui-master/backends/imgui_impl_glfw.cpp
        imgui-master/backends/imgui_impl_opengl3.cpp
)

file(GLOB IMPLOT_SOURCES
        implot-master/*.cpp
)

# Create the executable
add_executable(Project main.cpp ${IMGUI_SOURCES} ${IMPLOT_SOURCES})

# Link libraries
find_package(OpenGL REQUIRED)
find_package(glfw3 REQUIRED)

target_link_libraries(Project PRIVATE
        glfw
        OpenGL::GL
        dl
        pthread
)
