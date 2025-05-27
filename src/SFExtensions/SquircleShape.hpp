#ifndef SQUIRCLE_SHAPE_HPP
#define SQUIRCLE_SHAPE_HPP

#include <SFML/Graphics/Shape.hpp>
#include <cmath>

class SquircleShape : public sf::Shape {
    public:
        explicit SquircleShape (const sf::Vector2f& size = sf::Vector2f(0, 0), size_t pointPerCircle = 0, float radius = 0);
        ~SquircleShape ();
        void setSize (const sf::Vector2f& size);
        const sf::Vector2f& getSize () const;
        void setRadius (const float radius);
        float getRadius () const;
        virtual size_t getPointCount () const;
        virtual sf::Vector2f getPoint (size_t index) const;
    private:
        sf::Vector2f m_size;
        size_t m_pointCount;
        float m_radius;
};

#endif