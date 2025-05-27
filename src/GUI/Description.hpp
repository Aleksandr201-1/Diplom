#ifndef DESCRIPTION_HPP
#define DESCRIPTION_HPP

#include <GUI/GUIElement.hpp>
#include <GUI/TextBased.hpp>

class Description : public GUIElement, public TextBased {
    private:
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
    public:
        Description ();
        Description (const sf::String &str);
        ~Description ();
    private:
        sf::RectangleShape m_frame;
};

#endif