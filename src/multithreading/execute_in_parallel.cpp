#include "execute_in_parallel.hpp"

namespace multithreading {

void executeInParallel(std::vector<std::function<void()>> tasks) {
    int N_threads = std::thread::hardware_concurrency() - 1;
    // int N_threads = 1;
    boost::asio::thread_pool pool(N_threads);
    int N_tasks = tasks.size();
    for(int i = 0; i < N_tasks; ++i) {
        boost::asio::post(pool, tasks[i]);
    }
    pool.join();
}

}