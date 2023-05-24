#include <vector>
#include <string>
#include <fstream>

namespace utility {

template <typename T>
void writeCsv(std::string path, std::vector<T> data);

std::vector<double> readCsv(std::string path);

}