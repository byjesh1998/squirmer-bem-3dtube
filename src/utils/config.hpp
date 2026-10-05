// config.hpp -- tiny "key = value" input-file reader
//
// File format:
//   # comment
//   mode      = duct          # trailing comments allowed
//   a_over_R  = 0.3
//   dP_list   = 0.5 1 2 4     # lists are whitespace or comma separated
// Values given on the command line as key=value override the file.
#pragma once
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace bem {

class Config {
public:
    Config() = default;

    static Config load(const std::string& path) {
        std::ifstream in(path);
        if (!in) throw std::runtime_error("cannot open input file: " + path);
        Config c;
        std::string line;
        while (std::getline(in, line)) c.parse_line(line);
        return c;
    }
    void parse_line(std::string line) {
        auto h = line.find('#');
        if (h != std::string::npos) line = line.substr(0, h);
        auto eq = line.find('=');
        if (eq == std::string::npos) return;
        std::string k = trim(line.substr(0, eq)), v = trim(line.substr(eq + 1));
        if (!k.empty()) kv_[k] = v;
    }
    bool has(const std::string& k) const { return kv_.count(k) > 0; }
    std::string str(const std::string& k, const std::string& def = "") const {
        return has(k) ? kv_.at(k) : def;
    }
    double num(const std::string& k, double def) const { return has(k) ? std::stod(kv_.at(k)) : def; }
    double num(const std::string& k) const {
        if (!has(k)) throw std::runtime_error("missing required input: " + k);
        return std::stod(kv_.at(k));
    }
    int integer(const std::string& k, int def) const { return has(k) ? std::stoi(kv_.at(k)) : def; }
    std::vector<double> list(const std::string& k, std::vector<double> def = {}) const {
        if (!has(k)) return def;
        std::string s = kv_.at(k);
        for (char& c : s)
            if (c == ',') c = ' ';
        std::istringstream is(s);
        std::vector<double> v;
        double x;
        while (is >> x) v.push_back(x);
        return v;
    }
    const std::map<std::string, std::string>& all() const { return kv_; }


    
private:
    std::map<std::string, std::string> kv_;
    static std::string trim(const std::string& s) {
        auto b = s.find_first_not_of(" \t\r\n"), e = s.find_last_not_of(" \t\r\n");
        return b == std::string::npos ? "" : s.substr(b, e - b + 1);
    }
};

}  // namespace bem
