#ifndef ELLIPSE_SHAPE_HPP
#define ELLIPSE_SHAPE_HPP

#include <SFML/Graphics/Shape.hpp>
#include <cmath>

class EllipseShape : public sf::Shape {
    public:
        explicit EllipseShape (const sf::Vector2f& radius = sf::Vector2f(0, 0), size_t pointCount = 30);
        ~EllipseShape ();
        void setRadius (const sf::Vector2f& radius);
        const sf::Vector2f& getRadius () const;
        virtual size_t getPointCount () const;
        virtual sf::Vector2f getPoint (size_t index) const;
    private:
        sf::Vector2f m_radius;
        size_t m_pointCount;
};

#endif