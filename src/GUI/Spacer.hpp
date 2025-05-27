#ifndef SPACER_HPP
#define SPACER_HPP

#include <GUI/GUIElement.hpp>

class Spacer : public GUIElement {
    private:
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const;
    public:
        Spacer ();
        Spacer (float width, float height);
        Spacer (const sf::Vector2f &size);
        ~Spacer ();

        using GUIElement::setSize;
        void setSize (const sf::Vector2f &size) override;
        void saveToFile (std::ofstream &file) const override;
    private:
        static const std::string M_NAME;
        static const GuiFabric<Spacer> m_spacerFabric;
};

#endif