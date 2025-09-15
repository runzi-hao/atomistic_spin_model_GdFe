#ifndef IO_CSV_UTILS_H
#define IO_CSV_UTILS_H
#include <cctype>
#include <sstream>
#include <string>
#include <vector>

namespace io::csv {
    inline std::string trim(const std::string& s) {
        size_t a = 0;
        size_t b = s.size();
        while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
        while (a < b && std::isspace(static_cast<unsigned char>(s[b-1]))) --b;
        return s.substr(a, b-a);
    }

    inline std::vector<std::string> split_line_csv(const std::string& line) {
        std::vector<std::string> out;
        std::istringstream ss(line);
        std::string item;
        while (std::getline(ss, item, ','))
            out.push_back(trim(item));
        return out;
    }

} // namespace io::csv

#endif //IO_CSV_UTILS_H
