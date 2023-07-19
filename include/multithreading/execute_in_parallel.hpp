#include <boost/asio/thread_pool.hpp>
#include <boost/asio.hpp>
#include <functional>
#include <thread>
#include <vector>
#include <optional>

namespace multithreading {

void executeInParallel(std::vector<std::function<void()>> tasks);

template <typename ReturnType>
std::vector<ReturnType>
executeInParallel(std::vector<std::function<ReturnType()>> tasks) {
    int N_tasks = tasks.size();
    std::vector<ReturnType> output(N_tasks);
    std::vector<std::function<void()>> void_functions(N_tasks);
    for(int i = 0; i < N_tasks; ++i) {
        void_functions[i] = [&tasks, &output, i] () {
            output[i] = tasks[i]();
        };
    }
    executeInParallel(void_functions);
    return output;
}

template <typename ReturnType>
std::vector<std::optional<ReturnType>>
executeInParallelOpt(std::vector<std::function<std::optional<ReturnType>()>>
                     tasks) {
    int N_tasks = tasks.size();
    std::vector<std::optional<ReturnType>> output(N_tasks);
    std::vector<std::function<void()>> void_functions(N_tasks);
    for(int i = 0; i < N_tasks; ++i) {
        void_functions[i] = [&tasks, &output, i] () {
            output[i].emplace(tasks[i]().value());
        };
    }
    executeInParallel(void_functions);
    return output;
}



}