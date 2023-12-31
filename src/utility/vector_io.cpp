#include "vector_io.hpp"

namespace utility {

template <typename T>
void writeCsv(std::string path, std::vector<T> data) {
    std::ofstream fout(path);
    double N_lines = data.size();
    for(int i = 0; i < N_lines; ++i) {
        fout << data[i] << std::endl;
    }
    fout.close();
}

template void writeCsv(std::string path, std::vector<double> data);
// template void write_csv(std::string path, std::vector<int> data);


std::vector<double> readCsv(std::string path) {
    std::ifstream infile(path);
    std::string line;
    std::vector<double> output;
    while(std::getline(infile, line)) {
        output.push_back(std::stod(line));
    }
    infile.close();
    return output;
}

bool fileExists(std::string path) {
    // from https://www.geeksforgeeks.org/how-to-check-a-file-or-directory-exists-in-cpp/
    const char* file = path.c_str();
    struct stat sb;
    if(stat(file, &sb) == 0 && !(sb.st_mode & S_IFDIR)) {
        return 1;
    }
    else {
        return 0;
    }
}



}