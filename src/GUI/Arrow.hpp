#ifndef ARROW_HPP
#define ARROW_HPP

#include <GUI/GUIElement.hpp>
#include <initializer_list>
#include <SFExtensions/Line.hpp>
#include <memory>

class Arrow : public GUIElement {
    private:
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
    public:
        Arrow ();
        Arrow (const std::initializer_list<sf::Vector2f> &list);
        ~Arrow ();

        using GUIElement::setSize;
        void setSize (const sf::Vector2f &size) override;
        void loadFromFile (std::ifstream &file) override;
        void saveToFile (std::ofstream &file) const override;
    private:
        std::dynarray<sf::Vector2f> arr;
        float m_length;
        static const std::string M_NAME;
        static const GuiFabric<Arrow> m_arrowFabric;
};

#endif