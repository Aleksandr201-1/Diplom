#ifndef RHOMB_SHAPE_HPP
#define RHOMB_SHAPE_HPP

#include <SFML/Graphics/Shape.hpp>
#include <cmath>

class RhombShape : public sf::Shape {
    public:
        explicit RhombShape (const sf::Vector2f& size = sf::Vector2f(0, 0));
        ~RhombShape ();
        void setSize (const sf::Vector2f& size);
        const sf::Vector2f& getSize () const;
        virtual size_t getPointCount () const;
        virtual sf::Vector2f getPoint (size_t index) const;
    private:
        static const size_t M_POINT_COUNT = 4;
        sf::Vector2f m_size;
};

#endif