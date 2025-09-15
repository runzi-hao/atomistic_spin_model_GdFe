#include "io_temperature_csv.h"
#include "io_csv_utils.h"
#include <cctype>
#include <fstream>
#include <stdexcept>
using namespace io::csv;

void read_temperature_series_csv(const std::string& path,
    std::vector<double>& Te_kelvin)
{
    std::ifstream fin(path);
    if (!fin)
        throw std::runtime_error(
            "io_temperature_csv: Failed to open file: " + path);

    Te_kelvin.reserve(1024);
    std::string line;
    size_t line_no = 1;
    while (std::getline(fin, line)) {
        auto line_trimed = trim(line);
        if (line_trimed.empty() || line_trimed[0]=='#') continue;
        try {
            Te_kelvin.push_back(std::stod(line_trimed));
        }
        catch (...) {
            throw std::runtime_error(
                "Parse error at line " + std::to_string(line_no));
        }
        ++line_no;
    }
    if (Te_kelvin.empty())
        throw std::runtime_error("No data in file: " + path);
}

