#ifndef LINE_HPP
#define LINE_HPP

#include <SFML/Graphics.hpp>
#include <SFML/OpenGL.hpp>
#include <cmath>

class Line : public sf::Transformable, public sf::Drawable {
    public:
        enum class Type {
            POINT,
            LINE,
            LINE_WITH_POINT,
            DOTTED,
            DOTTED_WITH_POINT
        };
    private:
        void draw (sf::RenderTarget& target, sf::RenderStates states) const override;
    public:
        Line ();
        Line (const Line &line);
        Line (Type type, sf::Color color = sf::Color::Black, float thickness = 1.f, float interval = 0.f);
        ~Line ();

        void setInterval (float interval);

        Line &operator= (const Line &line);
    public:
        Type type;
        sf::Color color;
        float thickness, interval, intervalMiss;
        sf::Vector2f position[2];
    // private:
    //     float length = std::sqrt(xLen * xLen + yLen * yLen);
    //     intervalCount = (uint64_t)(length / interval);
};

class LineArray;

#endif