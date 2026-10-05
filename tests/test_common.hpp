// test_common.hpp -- tiny check/report helper shared by the tests
#pragma once
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>

namespace testing {

inline int& failures() { static int n = 0; return n; }
inline int& checks() { static int n = 0; return n; }

// value must satisfy |value - target| <= tol (absolute)
inline void check_close(const std::string& what, double value, double target, double tol) {
    ++checks();
    bool ok = std::fabs(value - target) <= tol;
    if (!ok) ++failures();
    std::printf("  [%s] %-58s value %.6g  target %.6g  (tol %.1e)\n", ok ? "PASS" : "FAIL",
                what.c_str(), value, target, tol);
}
inline void check_true(const std::string& what, bool ok, const std::string& detail = "") {
    ++checks();
    if (!ok) ++failures();
    std::printf("  [%s] %-58s %s\n", ok ? "PASS" : "FAIL", what.c_str(), detail.c_str());
}
inline bool has_flag(int argc, char** argv, const char* flag) {
    for (int i = 1; i < argc; ++i)
        if (!std::strcmp(argv[i], flag)) return true;
    return false;
}
inline int summary(const char* name) {
    std::printf("\n%s: %d/%d checks passed\n", name, checks() - failures(), checks());
    return failures() ? 1 : 0;
}

}  // namespace testing
