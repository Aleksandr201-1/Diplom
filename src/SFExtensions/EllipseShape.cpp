#include <SFExtensions/EllipseShape.hpp>

EllipseShape::EllipseShape (const sf::Vector2f& radius, size_t pointCount) : m_radius(radius), m_pointCount(pointCount) {
    update();
}

EllipseShape::~EllipseShape () {}

void EllipseShape::setRadius(const sf::Vector2f& radius) {
    m_radius = radius;
    update();
}

const sf::Vector2f& EllipseShape::getRadius() const {
    return m_radius;
}

size_t EllipseShape::getPointCount() const {
    return m_pointCount;
}

sf::Vector2f EllipseShape::getPoint(size_t index) const {
    static const float pi = std::acos(-1);

    float angle = index * 2 * pi / getPointCount() - pi / 2;
    float x = std::cos(angle) * m_radius.x;
    float y = std::sin(angle) * m_radius.y;

    return sf::Vector2f(m_radius.x + x, m_radius.y + y);
}