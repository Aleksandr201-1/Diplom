#include <General/Differentiation.hpp>

std::vector<std::pair<std::function<float128_t (float128_t)>, std::function<float128_t (float128_t)>>> func1;

int main () {
    std::cout.precision(20);
    std::cout.setf(std::ios::scientific);
    std::cout << "Its testing time\n";
    std::cout << "=====Correct result testing=====\n";
    func1 = {
        {
            [] (float128_t x) -> float128_t {
                return std::sin(x);
            },
            [] (float128_t x) -> float128_t {
                return std::cos(x);
            }
        },
        {
            [] (float128_t x) -> float128_t {
                return x*x*x;
            },
            [] (float128_t x) -> float128_t {
                return 3*x*x;
            }
        }
    };
    for (uint64_t i = 0; i < func1.size(); ++i) {
        float128_t diff = std::abs(derivative(func1[i].first, 5, 0.001, DiffConfig::POINTS2_ORDER1_WAY1) - func1[i].second(5));
        std::cout << "Test " << i + 1 << ": " << diff << "\n";
    }
    return 0;
}