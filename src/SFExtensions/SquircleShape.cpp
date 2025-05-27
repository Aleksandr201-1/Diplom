#include <SFExtensions/SquircleShape.hpp>

SquircleShape::SquircleShape (const sf::Vector2f& size, size_t pointPerCircle, float radius) : m_size(size), m_pointCount(pointPerCircle), m_radius(radius) {
    m_pointCount = m_pointCount * 4 + 8;
    float min = std::min(size.x, size.y);
    if (min / 2 < m_radius) {
        m_radius = min / 2;
    }
    update();
}

SquircleShape::~SquircleShape () {}

void SquircleShape::setSize(const sf::Vector2f& size) {
    m_size = size;
    update();
}

const sf::Vector2f& SquircleShape::getSize() const {
    return m_size;
}

void SquircleShape::setRadius (const float radius) {
    m_radius = radius;
    update();
}

float SquircleShape::getRadius () const {
    return m_radius;
}

size_t SquircleShape::getPointCount() const {
    return m_pointCount; // fixed, but could be an attribute of the class if needed
}

sf::Vector2f SquircleShape::getPoint(size_t index) const {
    static const float pi = std::acos(-1);
    size_t pointPerCircle = (getPointCount() - 8) / 4;
    if (index == 0) {
        return sf::Vector2f(m_radius, 0);
    } else if (index == 1) {
        return sf::Vector2f(m_size.x - m_radius, 0);
    } else if (index == pointPerCircle + 2) {
        return sf::Vector2f(m_size.x, m_radius);
    } else if (index == pointPerCircle + 3) {
        return sf::Vector2f(m_size.x, m_size.y - m_radius);
    } else if (index == pointPerCircle * 2 + 4) {
        return sf::Vector2f(m_size.x - m_radius, m_size.y);
    } else if (index == pointPerCircle * 2 + 5) {
        return sf::Vector2f(m_radius, m_size.y);
    } else if (index == pointPerCircle * 3 + 6) {
        return sf::Vector2f(0, m_size.y - m_radius);
    } else if (index == pointPerCircle * 3 + 7) {
        return sf::Vector2f(0, m_radius);
    }

    //return m_size / 2.f;

    float x = 0, y = 0;
    if (index < pointPerCircle + 2) {
        float angle = (index - 2) * 0.5 * pi / pointPerCircle;
        x = std::sin(angle) * m_radius - m_radius + m_size.x;
        y = -std::cos(angle) * m_radius + m_radius;
    } else if (index > pointPerCircle + 3 && index < pointPerCircle * 2 + 4) {
        float angle = (index - pointPerCircle - 4) * 0.5 * pi / pointPerCircle;
        x = std::cos(angle) * m_radius - m_radius + m_size.x;
        y = std::sin(angle) * m_radius + m_size.y - m_radius;
    } else if (index > pointPerCircle * 2 + 5 && index < pointPerCircle * 3 + 6) {
        float angle = (index - pointPerCircle * 2 - 6) * 0.5 * pi / pointPerCircle;
        x = -std::sin(angle) * m_radius + m_radius;
        y = std::cos(angle) * m_radius + m_size.y - m_radius;
    } else if (index > pointPerCircle * 3 + 7) {
        float angle = (index - pointPerCircle * 3 - 8) * 0.5 * pi / pointPerCircle;
        x = -std::cos(angle) * m_radius + m_radius;
        y = -std::sin(angle) * m_radius + m_radius;
    }

    return sf::Vector2f(x, y);
}