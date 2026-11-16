// Test oracle only. Prints libstdc++ random sequences for bit-exact replication checks.
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>
static bool first = true;
template <class V> static void vec(const std::string& k, const V& v) {
    std::cout << (first ? "" : ",\n") << "\"" << k << "\":[";
    first = false;
    std::size_t i = 0;
    for (auto x : v) std::cout << (i++ ? "," : "") << std::setprecision(17) << x;
    std::cout << "]";
}
int main() {
    std::cout << "{\n";
    std::mt19937_64 a(42);
    std::vector<std::string> raw;
    for (int i = 0; i < 700; ++i) { auto v = a(); if (i < 5 || i > 690) raw.push_back("\"" + std::to_string(v) + "\""); }
    vec("raw", raw);
    std::mt19937_64 b(7);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::vector<double> uu; for (int i = 0; i < 20; ++i) uu.push_back(u(b));
    vec("uniform", uu);
    std::uniform_int_distribution<std::size_t> pick(0, 6);
    std::vector<double> pk; for (int i = 0; i < 40; ++i) pk.push_back(pick(b));
    vec("pick", pk);
    std::uniform_int_distribution<std::size_t> pick2(3, 1000000007);
    std::vector<double> pk2; for (int i = 0; i < 10; ++i) pk2.push_back(pick2(b));
    vec("pick2", pk2);
    std::normal_distribution<double> nd;
    std::vector<double> nn; for (int i = 0; i < 21; ++i) nn.push_back(nd(b));
    vec("normal", nn);
    std::cout << "\n}\n";
}
