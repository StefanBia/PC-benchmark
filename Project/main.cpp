#include <cstdio>
#include "imgui-master/imgui.h"
#include "implot-master/implot.h"
#include "GLFW/glfw3.h"
#include "imgui-master/backends/imgui_impl_glfw.h"
#include "imgui-master/backends/imgui_impl_opengl3.h"
#include <cmath>
#include <array>
#include <memory>
#include <sstream>

std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

// GLFW error callback
void glfw_error_callback(int error, const char* description) {
    fprintf(stderr, "GLFW Error %d: %s\n", error, description);
}

int main() {
    // Set up GLFW
    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit()) return -1;

    // Create window with OpenGL context
    GLFWwindow* window = glfwCreateWindow(1280, 720, "CPU Model Info", NULL, NULL);
    if (!window) {
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); // Enable vsync

    // Set up Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;

    // Set up ImGui style
    ImGui::StyleColorsDark();

    // Initialize ImGui
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 130");

    std::string cpuModelLine = exec("lscpu | grep \"Model name\"");
    std::string cpuFreqMax = exec("lscpu | grep \"CPU max MHz\"");
    std::string cpuFreqMin = exec("lscpu | grep \"CPU min MHz\"");
    std::string cpuCores = exec("lscpu | grep \"Core(s)\"");
    std::string cpuThreads = exec("lscpu | grep \"Thread(s)\"");
    std::string totalRAM = exec("vmstat -s | grep \"total memory\"");

    // Main loop
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGui::Begin("CPU Information");
        ImGui::Text("%s", cpuModelLine.c_str());
        ImGui::Text("%s", cpuFreqMax.c_str());
        ImGui::Text("%s", cpuFreqMin.c_str());
        ImGui::Text("%s", cpuCores.c_str());
        ImGui::Text("%s", cpuThreads.c_str());
        ImGui::End();

        ImGui::Begin("Memory Information");
        ImGui::Text("%s", totalRAM.c_str());
        ImGui::End();

        // Render
        ImGui::Render();
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(0.45f, 0.55f, 0.60f, 1.00f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        // Swap buffers
        glfwSwapBuffers(window);
    }

    // Clean up
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}