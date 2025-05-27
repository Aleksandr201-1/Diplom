#ifndef POLYGON_SHAPE_HPP
#define POLYGON_SHAPE_HPP

#include <SFML/Graphics/Shape.hpp>
#include <cmath>

class PolygonShape : public sf::Shape {
    public:
        explicit PolygonShape (const sf::Vector2f& radius = sf::Vector2f(0, 0), size_t pointCount = 3);
        ~PolygonShape ();
        void setRadius (const sf::Vector2f& radius);
        const sf::Vector2f& getRadius () const;
        virtual size_t getPointCount () const;
        virtual sf::Vector2f getPoint (size_t index) const;
    private:
        sf::Vector2f m_radius;
        size_t m_pointCount;
};

#endif