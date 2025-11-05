// convert.cpp
// g++ -std=c++17 -o convert convert.cpp
#include <fstream>
#include <sstream>
#include <string>
#include <array>
#include <vector>
#include <cstdint>
#include <iostream>
#include <cctype>

static bool convert_dat_to_bin(const char* txtfile, const char* binfile)
{
    std::ifstream in(txtfile);
    if (!in.is_open()) {
        std::cerr << "Failed to open input: " << txtfile << "\n";
        return false;
    }

    // 1) scan to count rows and columns
    uint64_t nrows = 0;
    uint32_t ncols = 0;
    std::string line;
    while (std::getline(in, line)) {
        size_t p = 0;
        while (p < line.size() && std::isspace(static_cast<unsigned char>(line[p]))) ++p;
        if (p >= line.size() || line[p] == '#') continue;
        std::istringstream iss(line);
        std::vector<double> row;
        double v;
        while (iss >> v) row.push_back(v);
        if (row.empty()) continue;
        if (ncols == 0) ncols = static_cast<uint32_t>(row.size());
        else if (row.size() != ncols) {
            std::cerr << "Inconsistent columns at row " << (nrows + 1)
                      << ": expected " << ncols << ", got " << row.size() << "\n";
            return false;
        }
        ++nrows;
    }
    if (nrows == 0) {
        std::cerr << "No data in input: " << txtfile << "\n";
        return false;
    }
    in.clear();
    in.seekg(0, std::ios::beg);

    std::ofstream out(binfile, std::ios::binary);
    if (!out.is_open()) {
        std::cerr << "Failed to open output: " << binfile << "\n";
        return false;
    }

    uint64_t header_rows = 0;
    uint32_t header_cols = 0;
    out.write(reinterpret_cast<const char*>(&header_rows), sizeof(header_rows));
    out.write(reinterpret_cast<const char*>(&header_cols), sizeof(header_cols));

    const std::streamoff header_bytes = static_cast<std::streamoff>(sizeof(uint64_t) + sizeof(uint32_t));

    uint64_t row_idx = 0;
    while (std::getline(in, line)) {
        size_t p = 0;
        while (p < line.size() && std::isspace(static_cast<unsigned char>(line[p]))) ++p;
        if (p >= line.size() || line[p] == '#') continue;
        std::istringstream iss(line);
        std::vector<double> row;
        double v;
        while (iss >> v) row.push_back(v);
        if (row.empty()) continue;
        for (uint32_t j = 0; j < ncols; ++j) {
            std::streamoff pos = header_bytes + static_cast<std::streamoff>(sizeof(double) * (uint64_t(j) * nrows + row_idx));
            out.seekp(pos, std::ios::beg);
            out.write(reinterpret_cast<const char*>(&row[j]), sizeof(double));
        }
        ++row_idx;
    }

    out.seekp(0, std::ios::beg);
    out.write(reinterpret_cast<const char*>(&nrows), sizeof(nrows));
    out.write(reinterpret_cast<const char*>(&ncols), sizeof(ncols));
    out.close();
    in.close();
    return true;
}

double liner_inter(double x0, double y0, double x1, double y1, double x)
{
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}

// FastTable used for reading column-major binary files
const size_t COLUMN_COUNT = 4;
struct FastTable {
    std::vector<double> xs;
    std::vector<std::array<double,COLUMN_COUNT>> ys;

    bool load_bin(const char* binfile)
    {
        std::ifstream in(binfile, std::ios::binary);
        if (!in.is_open()) return false;
        uint64_t nrows = 0;
        uint32_t ncols = 0;
        in.read(reinterpret_cast<char*>(&nrows), sizeof(nrows));
        in.read(reinterpret_cast<char*>(&ncols), sizeof(ncols));
        if (nrows == 0 || ncols < 2) return false; 

        xs.resize(nrows);
        ys.assign(nrows, std::array<double, COLUMN_COUNT>{});

        in.read(reinterpret_cast<char*>(xs.data()), sizeof(double) * nrows);

        for (uint32_t col = 1; col <= COLUMN_COUNT; ++col) {
            std::vector<double> colbuf(nrows);
            if (col < ncols) {
                in.read(reinterpret_cast<char*>(colbuf.data()), sizeof(double) * nrows);
            } else {
                for (uint64_t i = 0; i < nrows; ++i) {
                    colbuf[i] = 0.0;
                }
            }
            for (uint64_t i = 0; i < nrows; ++i) {
                ys[i][col-1] = colbuf[i];
            }
        }

        // if columns more than COLUMN_COUNT + 1, skip remaining columns
        if (ncols > COLUMN_COUNT + 1) {
            uint64_t remaining = (uint64_t)(ncols - (COLUMN_COUNT + 1)) * nrows;
            in.seekg(static_cast<std::streamoff>(sizeof(double) * remaining), std::ios::cur);
        }

        in.close();
        return true;
    }
    
    double value(int col, double x) {
        if (col < 0 || col >= int(COLUMN_COUNT) || xs.empty()) return 0.0;
        if (x <= xs.front()) return ys.front()[col];
        if (x >= xs.back()) return ys.back()[col];

        // Binary search for the interval [xs[i], xs[i+1]] containing x
        size_t left = 0, right = xs.size() - 1;
        while (left + 1 < right) {
            size_t mid = left + (right - left) / 2;
            if (xs[mid] <= x)
                left = mid;
            else
                right = mid;
        }
        return liner_inter(xs[left], ys[left][col], xs[left+1], ys[left+1][col], x);
    }
};



static void debug_print_bin(const char* binfile)
{
    FastTable table;
    if (!table.load_bin(binfile)) {
        std::cerr << "Failed to load bin file: " << binfile << "\n";
        return;
    }
    std::cout.precision(16);
    for (size_t i = 0; i < table.xs.size(); ++i) {
        std::cout << table.xs[i];
        for (int j = 0; j < 5; ++j) {
            std::cout << " " << table.ys[i][j];
        }
        std::cout << "\n";
    }
    std::cout << "\nExample: Interpolated values at x = 1.5\n";
    double x = 1.5;
    for (int col = 0; col < int(COLUMN_COUNT); ++col) {
        double val = table.value(col, x);
        std::cout << "col " << col << ": " << val << "\n";
    }
}


int main(int argc, char** argv)
{
    if (argc == 2) {
        debug_print_bin(argv[1]);
        return 0;
    }
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " input.txt output.bin\n";
        return 1;
    }
    if (!convert_dat_to_bin(argv[1], argv[2])) {
        std::cerr << "Conversion failed\n";
        return 2;
    }
    std::cerr << "Converted " << argv[1] << " -> " << argv[2] << "\n";
    return 0;
}
