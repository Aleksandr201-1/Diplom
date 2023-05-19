#include <iostream>
#include <ctime>
#include <chrono>
//#include "../../src/General/FuncMaker.hpp"
#include <General/General.hpp>
#include <General/FuncMaker.hpp>

const uint64_t TEST_SIZE = 1'000'000;
const uint64_t TIME_TEST_NUM = -1;
//microseconds
using duration_t = std::chrono::milliseconds;

struct TestModule {
    FunctionalTree dFunc;
    std::function<float128_t (const std::vector<float128_t> &)> sFunc;
};

std::vector<TestModule> testInit () {
    return {
        {
            {"x + 2", std::vector<std::string>{"x"}},
            [] (const std::vector<float128_t> &X) {
                return X[0] + 2;
            }
        },
        {
            {"sin(cos(x)) + cos(cos(x))", std::vector<std::string>{"x"}},
            [] (const std::vector<float128_t> &X) {
                return std::sin(std::cos(X[0])) + std::cos(std::cos(X[0]));
            }
        },
        {
            {"sin(-x) + ln(e^10) + acos(-1)^x", std::vector<std::string>{"x"}},
            [] (const std::vector<float128_t> &X) {
                return std::sin(-X[0]) + std::log(std::exp(10)) + std::pow(std::acos(-1), X[0]);
            }
        },
        {
            {"(11 - x^3 * (3*x - 8)) / (12 * (x - 2)^2 * (x - 3))", std::vector<std::string>{"x"}},
            [] (const std::vector<float128_t> &X) {
                return (11 - X[0]*X[0]*X[0] * (3*X[0] - 8)) / (12 * (X[0] - 2)*(X[0]-2) * (X[0] - 3));
            }
        },
        {
            {"x * exp(1 / x)", std::vector<std::string>{"x"}},
            [] (const std::vector<float128_t> &X) {
                return X[0] * exp(1.0 / X[0]);
            }
        },
        {
            {"exp(x*sin(ln(x)) + x)", std::vector<std::string>{"x"}},
            [] (const std::vector<float128_t> &X) {
                return std::exp(X[0]*std::sin(std::log(X[0])) + X[0]);
            }
        },
        {
            {"(1/4)*x^3 - 1/x", std::vector<std::string>{"x"}},
            [] (const std::vector<float128_t> &X) {
                return (1.0/4.0)*X[0]*X[0]*X[0] - 1.0/X[0];
            }
        },
        {
            {"1 + ln(abs(x))", std::vector<std::string>{"x"}},
            [] (const std::vector<float128_t> &X) {
                return 1 + std::log(std::abs(X[0]));
            }
        },
        {
            {"cos(sqrt(4*x)) + sin(sqrt(4*x))", std::vector<std::string>{"x"}},
            [] (const std::vector<float128_t> &X) {
                return std::cos(std::sqrt(4*X[0])) + std::sin(std::sqrt(4*X[0]));
            }
        },
        {
            {"abs(x)^(3/2)", std::vector<std::string>{"x"}},
            [] (const std::vector<float128_t> &X) {
                return std::pow(std::abs(X[0]), (3.0/2.0));
            }
        }
    };
}

int main () {
    std::cout.precision(20);
    std::cout.setf(std::ios::scientific);
    std::srand((uint64_t)std::time(nullptr));
    std::cout << "Its testing time\n";
    std::cout << "=====Correct result testing=====\n";

    auto test = testInit();
    float128_t diff, ans1, ans2;
    float128_t random;
    uint64_t dTime, sTime;
    uint64_t left, right;
    uint64_t errorSpot = -1;
    std::chrono::time_point <std::chrono::system_clock> startT, endT;

    for (uint64_t i = 0; i < test.size(); ++i) {
        diff = 0;
        errorSpot = -1;
        for (uint64_t j = 0; j < TEST_SIZE; ++j) {
            random = j * 1.0;
            ans1 = test[i].dFunc({random});
            ans2 = test[i].sFunc({random});
            if (std::isnan(ans1) && std::isnan(ans2)) {
                continue;
            }
            if (std::abs(ans1) == INFINITY && std::abs(ans2) == INFINITY && ans1 * ans2 > 0) {
                continue;
            }
            diff += std::abs(ans1 - ans2);
            if (diff != 0.0) {
                errorSpot = j;
                break;
            }
        }
        std::cout << "Test " << i + 1 << ": ";

        if (errorSpot != -1) {
            random = errorSpot * 1.0;
            std::cout << "Not OK\n";
            std::cout << "\tTotal diff: " << diff << "\n";
            std::cout << "\tAt the point: " << random << "\n";
            std::cout << "\tResult at " << random << "\n";
            std::cout << "\t\tTree(" << random << ") = " << test[i].dFunc({random}) << "\n";
            std::cout << "\t\tFunc(" << random << ") = " << test[i].sFunc({random}) << "\n";
            std::cout << "\tFunction: " << test[i].dFunc << "\n";
            random = (rand() % 100) * 1.0;
            std::cout << "\tChoosing random point: " << random << "\n";
            std::cout << "\tResult at " << random << "\n";
            std::cout << "\t\tTree(" << random << ") = " << test[i].dFunc({random}) << "\n";
            std::cout << "\t\tFunc(" << random << ") = " << test[i].sFunc({random}) << "\n";
        } else {
            std::cout << "OK\n";
        }
    }
    std::cout << "=====Time testing=====\n";
    if (TIME_TEST_NUM == -1) {
        std::cout << "Test num for time testing not set. Testing for all\n\n";
        left = 0;
        right = test.size();
    } else {
        left = TIME_TEST_NUM - 1;
        right = TIME_TEST_NUM;
    }
    for (uint64_t i = left; i < right; ++i) {
        std::cout << "Using test №" << i + 1 << "\n";
        dTime = sTime = 0;

        startT = std::chrono::system_clock::now();
        for (uint64_t j = 0; j < TEST_SIZE; ++j) {
            random = test[i].dFunc({j * 1.0});
        }
        endT = std::chrono::system_clock::now();
        dTime += std::chrono::duration_cast<duration_t>(endT - startT).count();

        startT = std::chrono::system_clock::now();
        for (uint64_t j = 0; j < TEST_SIZE; ++j) {
            random = test[i].sFunc({j * 1.0});
        }
        endT = std::chrono::system_clock::now();
        sTime += std::chrono::duration_cast<duration_t>(endT - startT).count();

        std::cout << "Dynamic func time: " << dTime << "mcs\n";
        std::cout << "Static func time:  " << sTime << "mcs\n\n";
    }
    std::cout << "Looks like we tested all other the place\n";
    return 0;
}