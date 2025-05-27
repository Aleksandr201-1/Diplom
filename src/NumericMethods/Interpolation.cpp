#include <NumericMethods/Interpolation.hpp>

double LinearInterpolation (const std::vector<std::pair<double, double>> &points, double x) {
    double ans = 0;
    uint64_t i = 0;
    while (x > points[i + 1].first) {
        ++i;
        if (i + 1 == points.size()) {
            --i;
            break;
        }
    }
    double coeff = (points[i + 1].second - points[i].second) / (points[i + 1].first - points[i].first);
    return (points[i].second - points[i].first * coeff) + coeff * x;
}

double LinearInterpolation (const std::vector<sf::Vector2f> &points, double x) {
    double ans = 0;
    uint64_t i = 0;
    while (x > points[i + 1].x) {
        ++i;
        if (i + 1 == points.size()) {
            --i;
            break;
        }
    }
    double coeff = (points[i + 1].y - points[i].y) / (points[i + 1].x - points[i].x);
    return (points[i].y - points[i].x * coeff) + coeff * x;
}