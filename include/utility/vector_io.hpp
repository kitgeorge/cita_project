#include <vector>
#include <string>
#include <fstream>
#include <sys/stat.h>

namespace utility {

template <typename T>
void writeCsv(std::string path, std::vector<T> data);

std::vector<double> readCsv(std::string path);

bool fileExists(std::string path);

}