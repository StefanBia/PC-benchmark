#include <cstdio>
#include "imgui-master/imgui.h"
#include "implot-master/implot.h"
#include "GLFW/glfw3.h"
#include "imgui-master/backends/imgui_impl_glfw.h"
#include "imgui-master/backends/imgui_impl_opengl3.h"
#include "implot-master/implot.h"
#include "implot-master/implot_internal.h"
#include <cmath>
#include <array>
#include <memory>
#include <sstream>
#include <x86intrin.h>
#include <cstdio>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>
#include <ctime>
#include <unistd.h>
#include <fcntl.h>
#include <dirent.h>
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <vector>


std::ostringstream results_read;
std::ostringstream results_write;
std::ostringstream results_cpu;
GLFWwindow* window;
float progress_read, progress_write, progress_cpu;
bool is_read_ready = false, show_progress_read = false, show_progress_write = false, show_progress_cpu;
bool is_write_ready = false, is_cpu_ready = false;
volatile bool running_thread = false;

///////////////////////////////////////////////////////////////////////////////////////////////////////

float bar_data[] = {32.11f,21.38, 0.0f, 13.97f};   // Bar heights as floats
float rw_data[] = {60.34, 35.23, 0.0, 108.73};
float ram_data[] = {760.09, 1426.71, 0.0, 2004.67};
const char* bar_labels[] = {"Intel Core i3-1115G4/HDD","Intel Core i5-1135G7/HDD", "Your PC", "AMD Ryzen 7 4700U/SSD"}; // Labels for bars
float x_positions[] = {0.0f, 1.0f, 2.0f, 3.0};

float max_bar_height = 0.0f;

float rd_data = 0.0;
float wr_data = 0.0;

float ram_rd_data = 0.0;
float ram_wr_data = 0.0;

//////////////////////////////////////////////////////////////////////////////////////////////////////

volatile float total_exec_time = 0.0;

inline uint64_t rdtscp() {
    unsigned int aux;
    return __rdtscp(&aux);
}

void test_arithmetic_and_logical_operations() {

    const int runs = 1;

#define MEASURE_OPERATION(name, asm_code)                             \
    {                                                                 \
        volatile int result = 0;                                      \
        uint64_t total_cycles = 0;                                    \
        for (int run = 0; run < runs; ++run) {                        \
            uint64_t start_cycles = rdtscp();                         \
            asm volatile (asm_code                                    \
                : "=a" (result)                                       \
                : "a" (result), "b" (3), "c" (2)                      \
                : "cc");                                              \
            uint64_t end_cycles = rdtscp();                           \
            total_cycles += (end_cycles - start_cycles);              \
        }                                                             \
        results_cpu << name << " cycles (avg): "                      \
                    << (total_cycles / runs) << "\n";            \
    }

#define MEASURE_FP_DIVISION(name)                                     \
    {                                                                 \
        volatile double result = 0.0;                                 \
        double dividend = static_cast<double>(30.56);                 \
        double divisor = static_cast<double>(9.354);                  \
        uint64_t total_cycles = 0;                                    \
        for (int run = 0; run < runs; ++run) {                        \
            uint64_t start_cycles = rdtscp();                         \
            asm volatile (                                            \
                "fldl %1;"                                            \
                "fldl %2;"                                            \
                "fdivp;"                                              \
                "fstpl %0;"                                           \
                : "=m" (result)                                       \
                : "m" (dividend), "m" (divisor)                       \
                : "st");                                              \
            uint64_t end_cycles = rdtscp();                           \
            total_cycles += (end_cycles - start_cycles);              \
        }                                                             \
        results_cpu << name << " cycles (avg): "                      \
                    << (total_cycles / runs) << "\n";            \
    }

    // Measure each operation
    MEASURE_FP_DIVISION("Floating-point Division");
    MEASURE_OPERATION("Addition", "add %%ebx, %%eax;")
    MEASURE_OPERATION("Subtraction", "sub %%ebx, %%eax;")
    MEASURE_OPERATION("Multiplication", "imul %%ebx, %%eax;")

    MEASURE_OPERATION("Bitwise AND", "and %%ebx, %%eax;")
    MEASURE_OPERATION("Bitwise OR", "or %%ebx, %%eax;")
    MEASURE_OPERATION("Bitwise XOR", "xor %%ebx, %%eax;")
    MEASURE_OPERATION("Left Shift", "shl $1, %%eax;")
    MEASURE_OPERATION("Right Shift", "shr $1, %%eax;")
}





struct Point {
    int x;
    int y;
};

int a = 1;
int b = 1;

int mod_inverse(int a, int mod) {
    int t = 0, new_t = 1;
    int r = mod, new_r = a % mod;

    while (new_r != 0) {
        int quotient = r / new_r;
        t -= quotient * new_t;
        std::swap(t, new_t);
        r -= quotient * new_r;
        std::swap(r, new_r);
    }

    if (r > 1) return -1;
    if (t < 0) t += mod;

    return t;
}

Point elliptic_add(const Point& P, const Point& Q, int N) {
    if (P.x == Q.x && P.y == Q.y) {
        int lambda = ((3 * P.x * P.x + a) * mod_inverse(2 * P.y, N)) % N;
        int x = (lambda * lambda - 2 * P.x) % N;
        int y = (lambda * (P.x - x) - P.y) % N;
        return Point{x, y};
    }

    int lambda = ((Q.y - P.y) * mod_inverse(Q.x - P.x, N)) % N;
    int x = (lambda * lambda - P.x - Q.x) % N;
    int y = (lambda * (P.x - x) - P.y) % N;

    return Point{x, y};
}

Point elliptic_multiply(int k, const Point& P, int N) {
    Point result = {0, 0};
    Point current = P;

    while (k > 0) {
        if (k % 2 == 1) {
            result = elliptic_add(result, current, N);
        }
        current = elliptic_add(current, current, N);
        k /= 2;
    }

    return result;
}

int ecm_factor(int N, int max_attempts = 100) {
    for (int attempt = 0; attempt < max_attempts; ++attempt) {
        a = rand() % N;
        b = rand() % N;

        Point P = {rand() % N, rand() % N};

        int k = 2;
        while (true) {
            Point Q = elliptic_multiply(k, P, N);
            if (Q.x == 0 && Q.y == 0) {
                break;
            }
            k += 1;

            int factor = std::__gcd(Q.x, N);
            if (factor > 1 && factor < N) {
                break;
            }
        }
    }
    return -1;
}

bool ifnotPrime(int prime[], int x)
{
    return (prime[x/64] & (1 << ((x >> 1) & 31)));
}

bool makeComposite(int prime[], int x)
{
    prime[x/64] |= (1 << ((x >> 1) & 31));
}


void bitWiseSieve(int n)
{
    int prime[n/64];
    memset(prime, 0, sizeof(prime));

    for (int i = 3; i * i <= n; i += 2) {
        if (!ifnotPrime(prime, i))
            for (int j = i * i, k = i << 1; j < n; j += k)
                makeComposite(prime, j);
    }

    for (int i = 3; i <= n; i += 2)
        if (!ifnotPrime(prime, i)){
            1;
        }

}
void pi_decimals(int digits){
    double pi;
        for (int k = 0; k < digits; k++) {
        double term1 = 1.0 / pow(16, k) * (4.0 / (8 * k + 1));
        double term2 = 1.0 / pow(16, k) * (2.0 / (8 * k + 4));
        double term3 = 1.0 / pow(16, k) * (1.0 / (8 * k + 5));
        double term4 = 1.0 / pow(16, k) * (1.0 / (8 * k + 6));

        pi += term1 - term2 - term3 - term4;
    }
}

// Function to measure execution time of arithmetic and logical operations
void test_cpu() {
    int digits = 100000000;
    int limit = 100000000;


    auto start_time = std::chrono::high_resolution_clock::now();

    pi_decimals(1000000);
    progress_cpu += 1.0/7;
    pi_decimals(10000000);
    progress_cpu += 1.0/7;
    pi_decimals(100000000);
    progress_cpu += 1.0/7;

    auto end_time = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end_time - start_time;
    results_cpu << "Elapsed time: " << elapsed.count() << "s for calculating decimals of PI" << std::endl;
    total_exec_time = elapsed.count();
    start_time = std::chrono::high_resolution_clock::now();

    bitWiseSieve(1000000);
    progress_cpu += 1.0/7;
    bitWiseSieve(10000000);
    progress_cpu += 1.0/7;
    bitWiseSieve(100000000);
    progress_cpu += 1.0/7;

    end_time = std::chrono::high_resolution_clock::now();

    elapsed = end_time - start_time;
    results_cpu << "Elapsed time: " << elapsed.count() << " s for Bitwise Eratosthenes Sieve" << std::endl;
    total_exec_time += elapsed.count();
    int N = 1999999999;
    start_time = std::chrono::high_resolution_clock::now();

    ecm_factor(N, 100000);

    end_time = std::chrono::high_resolution_clock::now();
    progress_cpu += 1.0/7;
    elapsed = end_time - start_time;
    results_cpu << "Elapsed time: " << elapsed.count() << " s for Lenstra Elliptic-Curve Factorization" << std::endl;
    total_exec_time += elapsed.count();

    test_arithmetic_and_logical_operations();
    is_cpu_ready = true;
    show_progress_cpu = false;
    progress_cpu = 0.0;
    running_thread = false;

    bar_data[2] = total_exec_time;

    for (int i = 0; i < 3; i++) {
        if (bar_data[i] > max_bar_height) {
            max_bar_height = bar_data[i];
        }
    }

}





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

void measure_write_speeds(){
    int fd = open("../test-files/write.txt", O_WRONLY | O_CREAT | O_DIRECT, 0644);
    if (fd < 0) {
        perror("error opening file");
        return;
    }

    const size_t BUFFER_SIZE = 4096;//4KB

    double mbps1 = 0.0, mbps5 = 0.0, mbps100 = 0.0;
    std::chrono::duration<double> elapsed;
    void* buffer;
    if (posix_memalign(&buffer, 4096, BUFFER_SIZE) != 0) {
        perror("memory allocation failed");
        close(fd);
        return;
    }

    int nr_runs = 5;

    for(int i = 0; i < nr_runs; i++) {

        auto start_time = std::chrono::high_resolution_clock::now();

        for (int i = 0; i <= 1250; i++) {///5MB
            write(fd, buffer, BUFFER_SIZE);
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        elapsed = end_time - start_time;
        mbps5 += (1250 * 1024 / 1048576.0) / elapsed.count();

        lseek(fd, 0, SEEK_SET);

        start_time = std::chrono::high_resolution_clock::now();

        for (int i = 0; i <= 25000; i++) {///100MB
            write(fd, buffer, BUFFER_SIZE);
        }

        end_time = std::chrono::high_resolution_clock::now();
        elapsed = end_time - start_time;
        mbps100 += (25000 * 1024 / 1048576.0) / elapsed.count();

        lseek(fd, 0, SEEK_SET);

        start_time = std::chrono::high_resolution_clock::now();

        for (int i = 0; i <= 250000; i++) {///1GB
            write(fd, buffer, BUFFER_SIZE);
        }

        end_time = std::chrono::high_resolution_clock::now();
        elapsed = end_time - start_time;
        mbps1 += (250000 * 1024 / 1048576.0) / elapsed.count();

        lseek(fd, 0, SEEK_SET);
        progress_write += 1.0 / nr_runs;

    }

    unlink("../test-files/write.txt");
    free(buffer);

    mbps5 /= nr_runs;
    mbps1 /= nr_runs;
    mbps100 /= nr_runs;

    wr_data = (mbps1 + mbps5 + mbps100) / 3;
    rw_data[2] = (rd_data + wr_data) / 2;

    results_write << mbps5 << "MB/s for 5MB file size" << std::endl;
    results_write << mbps100 << "MB/s for 100MB file size" << std::endl;
    results_write << mbps1 << "MB/s for 1GB file size" << std::endl;

    std::vector<char> memory((1<<26), 0); // 64MB

    auto start = std::chrono::high_resolution_clock::now();

    volatile char sink = 'a';
    for (size_t i = 0; i < memory.size(); i++) {
        memory[i] = sink;
    }

    auto end = std::chrono::high_resolution_clock::now();

    auto start1 = std::chrono::high_resolution_clock::now();

    int size = 1<<26;

    for (size_t i = 0; i < size; i++) {
        sink = 1;
    }

    auto end1 = std::chrono::high_resolution_clock::now();

    double duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    double duration_for = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1).count();
    duration_ms -= duration_for;
    double throughput = (memory.size() / (1024.0 * 1024.0)) / (duration_ms / 1000.0);

    ram_wr_data = throughput;
    ram_data[2] = (ram_wr_data + ram_rd_data)/2;
    results_write << "CPU-to-RAM transfer speed:" << throughput << " MB/s" << std::endl;

    show_progress_write = false;
    is_write_ready = true;
    progress_write = 0.0;
    running_thread = false;

}

std::string measure_data_transfer_speeds(){

    int fd = open("../../../test-files/file1mb.txt", O_RDONLY | O_DIRECT);
    if (fd < 0) {
        perror("error opening file");
        return nullptr;
    }
    int nr_runs = 5;
    off_t totalFileSize;
    struct stat fileStat;
    //////////////////////////////////////////////////////////
    ///Test Run
    const size_t BUFFER_SIZE = 4096;

    void* buffer;
    if (posix_memalign(&buffer, 4096, BUFFER_SIZE) != 0) {
        perror("memory allocation failed");
        close(fd);
        return nullptr;
    }

    auto start_time = std::chrono::high_resolution_clock::now();
    ssize_t bytesRead = 0, totalBytes = 0;

    while ((bytesRead = read(fd, buffer, BUFFER_SIZE)) > 0) {
        totalBytes += bytesRead;
    }

    auto end_time = std::chrono::high_resolution_clock::now();

    double mbps;
    std::chrono::duration<double> elapsed = end_time - start_time;

    mbps = (totalBytes / 1048576.0) / elapsed.count();
//    results << elapsed.count() << "s, " << mbps << " MB/s Read Speed for " << (totalBytes / 1024) << " KB file" << std::endl;
    free(buffer);

    double mbps1 = 0.0, mbps5 = 0.0, mbps100 = 0.0;


    for(int i = 0; i < nr_runs; i++) {

        ///////////////////////////////////////////////////////////////////////////////////////////
        ///5MB file

        fd = open("../../../test-files/file5mb.txt", O_RDONLY | O_DIRECT);

        if (fd < 0) {
            perror("Error opening file");
            return nullptr;
        }

        if (posix_memalign(&buffer, 4096, BUFFER_SIZE) != 0) {
            perror("memory allocation failed");
            close(fd);
            return nullptr;
        }

        start_time = std::chrono::high_resolution_clock::now();
        bytesRead = 0, totalBytes = 0;

        if (fstat(fd, &fileStat) < 0) {
            perror("Error getting file size");
            close(fd);
            return nullptr;
        }
        totalFileSize = fileStat.st_size;

        while ((bytesRead = read(fd, buffer, BUFFER_SIZE)) > 0) {
            totalBytes += bytesRead;
            float progress = (float)totalBytes / totalFileSize;

        }

        end_time = std::chrono::high_resolution_clock::now();


        elapsed = end_time - start_time;

        mbps5 += (totalBytes / 1048576.0) / elapsed.count();
//    results << elapsed.count() << "s, " << mbps << " MB/s Read Speed for " << (totalBytes / 1024) << " KB file" << std::endl;
        free(buffer);

        ///////////////////////////////////////////////////////////////////////////////////////////
        ///100MB file

        fd = open("../../../test-files/100MB.zip", O_RDONLY | O_DIRECT);

        if (fd < 0) {
            perror("Error opening file");
            return nullptr;
        }

        if (posix_memalign(&buffer, 4096, BUFFER_SIZE) != 0) {
            perror("memory allocation failed");
            close(fd);
            return nullptr;
        }

        start_time = std::chrono::high_resolution_clock::now();
        bytesRead = 0, totalBytes = 0;

        while ((bytesRead = read(fd, buffer, BUFFER_SIZE)) > 0) {
            totalBytes += bytesRead;
        }

        end_time = std::chrono::high_resolution_clock::now();


        elapsed = end_time - start_time;

        mbps100 += (totalBytes / 1048576.0) / elapsed.count();
//    results << elapsed.count() << "s, " << mbps << " MB/s Read Speed for " << (totalBytes / 1024) << " KB file" << std::endl;
        free(buffer);

        ///////////////////////////////////////////////////////////////////////////////////////////
        ///1GB file

        fd = open("../../../test-files/1GB.zip", O_RDONLY | O_DIRECT);

        if (fd < 0) {
            perror("Error opening file");
            return nullptr;
        }

        if (posix_memalign(&buffer, 4096, BUFFER_SIZE) != 0) {
            perror("memory allocation failed");
            close(fd);
            return nullptr;
        }

        start_time = std::chrono::high_resolution_clock::now();
        bytesRead = 0, totalBytes = 0;

        while ((bytesRead = read(fd, buffer, BUFFER_SIZE)) > 0) {
            totalBytes += bytesRead;
        }

        end_time = std::chrono::high_resolution_clock::now();


        elapsed = end_time - start_time;

        mbps1 += (totalBytes / 1048576.0) / elapsed.count();
//    results << elapsed.count() << "s, " << mbps << " MB/s Read Speed for " << (totalBytes / 1024) << " KB file" << std::endl;
        free(buffer);

        progress_read += 1.0 / nr_runs;

    }

    mbps1 /= nr_runs;
    mbps5 /= nr_runs;
    mbps100 /= nr_runs;
    rd_data = (mbps1 + mbps5 + mbps100) / 3;
    rw_data[2] = (rd_data + wr_data) / 2;

    results_read << mbps5 << " MB/s Read Speed for 5 MB file" << std::endl;
    results_read << mbps100 << " MB/s Read Speed for 100 MB file" << std::endl;
    results_read << mbps1 << " MB/s Read Speed for 1 GB file" << std::endl;

    std::vector<char> memory((1<<26), 0); // 64MB

    auto start = std::chrono::high_resolution_clock::now();

    volatile char sink;
    for (size_t i = 0; i < memory.size(); i++) {
        sink = memory[i];
    }

    auto end = std::chrono::high_resolution_clock::now();

    auto start1 = std::chrono::high_resolution_clock::now();

    int size = 1<<26;

    for (size_t i = 0; i < size; i++) {
        sink = 1;
    }

    auto end1 = std::chrono::high_resolution_clock::now();

    double duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    double duration_for = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1).count();
    duration_ms -= duration_for;
    double throughput = (memory.size() / (1024.0 * 1024.0)) / (duration_ms / 1000.0);


    ram_rd_data = throughput;
    ram_data[2] = (ram_wr_data + ram_rd_data)/2;
    results_read << "RAM-to-CPU transfer speed:" << throughput << " MB/s" << std::endl;

    std::string finalString = results_read.str();
    is_read_ready = true;
    show_progress_read = false;
    progress_read = 0.0;
    running_thread = false;
    return finalString;
}

void startDataThread(){
    pthread_t tid;
    if (pthread_create(&tid, nullptr, [](void*) -> void* {
        measure_data_transfer_speeds();
    }, nullptr) != 0) {
        perror("Error creating thread");
        exit(1);
    }


}

void startDataWriteThread(){
    pthread_t tid;
    if (pthread_create(&tid, nullptr, [](void*) -> void* {
        measure_write_speeds();
    }, nullptr) != 0) {
        perror("Error creating thread");
        exit(1);
    }

}

void startCpuTestThread(){
    pthread_t tid;
    if (pthread_create(&tid, nullptr, [](void*) -> void* {
        test_cpu();
    }, nullptr) != 0) {
        perror("Error creating thread");
        exit(1);
    }

}

// GLFW error callback
void glfw_error_callback(int error, const char* description) {
    fprintf(stderr, "GLFW Error %d: %s\n", error, description);
}



int main() {

    for (int i = 0; i < 3; i++) {
        if (bar_data[i] > max_bar_height) {
            max_bar_height = bar_data[i];
        }
    }

    // Set up GLFW
    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit()) return -1;

    // Create window with OpenGL context
    window = glfwCreateWindow(1920, 1080, "CPU Model Info", NULL, NULL);
    if (!window) {
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); // Enable vsync


    // Set up Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
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

        ImGui::Begin("TestCPU");
        if (ImGui::Button("Test CPU")&& !running_thread) {
            running_thread = true;
            is_cpu_ready = false;
            show_progress_cpu = true;
            startCpuTestThread();
        }
        if (ImGui::Button("Test CPU One Thread")&& !running_thread) {
            running_thread = true;
            is_cpu_ready = false;
            show_progress_cpu = true;
            test_cpu();
        }
        ImGui::Text("CPU speeds: ");

        if(is_cpu_ready){
            ImGui::Text("%s", results_cpu.str().c_str());
        }


        ImGui::End();

        ImGui::Begin("Read Speeds");

        if (ImGui::Button("Measure Reading Speed") && !running_thread) {
            running_thread = true;
            results_read.clear();
            startDataThread();
            show_progress_read = true;
        }
        if (ImGui::Button("Measure Reading Speed One Thread") && !running_thread) {
            running_thread = true;
            results_read.clear();
            measure_data_transfer_speeds();
            show_progress_read = true;
        }

        ImGui::Text("Read transfer speeds: ");

        if(is_read_ready)
            ImGui::Text("%s", results_read.str().c_str());

        ImGui::End();

        if(show_progress_read){
            ImGui::Begin("Read Progress");
            ImGui::Text("Testing read speed...");
            ImGui::ProgressBar(progress_read, ImVec2(0.0f, 0.0f), "Reading...");
            ImGui::End();
        }

        ImGui::Begin("Write Speeds");

        if (ImGui::Button("Measure Writing Speed") && !running_thread) {
            running_thread = true;
            results_write.clear();
            show_progress_write = true;
            startDataWriteThread();
        }
        if (ImGui::Button("Measure Writing Speed One Thread") && !running_thread) {
            running_thread = true;
            results_write.clear();
            show_progress_write = true;
            measure_write_speeds();
        }

        ImGui::Text("Write transfer speeds: ");

        if(is_write_ready)
            ImGui::Text("%s", results_write.str().c_str());

        ImGui::End();

        if(show_progress_write){
            ImGui::Begin("Write Progress");
            ImGui::Text("Testing write speed...");
            ImGui::ProgressBar(progress_write, ImVec2(0.0f, 0.0f), "Writing...");
            ImGui::End();
        }

        if(show_progress_cpu){
            ImGui::Begin("CPU Test Progress");
            ImGui::Text("Testing CPU...");
            ImGui::ProgressBar(progress_cpu, ImVec2(0.0f, 0.0f), "Testing...");
            ImGui::End();
        }

//        ImGui::Begin("graph");
//        if(ImPlot::BeginPlot("MyPlot")){
//            ImPlot::PlotBars("Bars", bar_data, 3);
//            ImPlot::EndPlot();
//        }
//        ImGui::End();
        ImGui::Begin("CpuGraph");
        ImPlot::SetNextAxesLimits(-2.5, 4.5, -0.5, max_bar_height + 1.0, ImPlotCond_Always);
        if (ImPlot::BeginPlot("CPU Execution Time")) {
            // Plot each bar with a unique label
            for (int i = 0; i < 4; i++) {
                // Use float pointer for the single bar
                ImPlot::PlotBars(bar_labels[i], &bar_data[i], 1, 0.5, x_positions[i]);
            }
            ImPlot::EndPlot();
        }
        ImGui::End();

        ImGui::Begin("RWGraph");
        ImPlot::SetNextAxesLimits(-2.5, 4.5, -0.5, 110.0, ImPlotCond_Always);
        if (ImPlot::BeginPlot("Read/Write Speeds")) {
            // Plot each bar with a unique label
            for (int i = 0; i < 4; i++) {
                // Use float pointer for the single bar
                ImPlot::PlotBars(bar_labels[i], &rw_data[i], 1, 0.5, x_positions[i]);
            }
            ImPlot::EndPlot();
        }
        ImGui::End();

        ImGui::Begin("RAMGraph");
        ImPlot::SetNextAxesLimits(-2.5, 4.5, -0.5, 2025.0, ImPlotCond_Always);
        if (ImPlot::BeginPlot("RAM Speeds")) {
            // Plot each bar with a unique label
            for (int i = 0; i < 4; i++) {
                // Use float pointer for the single bar
                ImPlot::PlotBars(bar_labels[i], &ram_data[i], 1, 0.5, x_positions[i]);
            }
            ImPlot::EndPlot();
        }
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
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}