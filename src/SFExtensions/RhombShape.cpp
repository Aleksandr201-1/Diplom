#include <SFExtensions/RhombShape.hpp>

RhombShape::RhombShape (const sf::Vector2f& size) : m_size(size) {}

RhombShape::~RhombShape () {}

void RhombShape::setSize (const sf::Vector2f& size) {
    m_size = size;
}

const sf::Vector2f& RhombShape::getSize () const {
    return m_size;
}

size_t RhombShape::getPointCount () const {
    return M_POINT_COUNT;
}

sf::Vector2f RhombShape::getPoint (size_t index) const {
    sf::Vector2f point;
    int signedId = static_cast<int>(index);
    point.x = m_size.x / 2 * (2 - std::abs(signedId - 2)); //0 1 2 1
    point.y = m_size.y / 2 * (std::abs(signedId - 1)); //1 0 1 2
    return point;
}