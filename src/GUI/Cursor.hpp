#ifndef CURSOR_HPP
#define CURSOR_HPP

#include <GUI/Interactable.hpp>

class Cursor : public Interactable {
    private:
        void updateContent () override;
        void interact (const sf::Event &event, const sf::Vector2f &mousePos) override;
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
    public:
        Cursor ();
        ~Cursor ();

        using GUIElement::setSize;
        void setSize (const sf::Vector2f &size) override;
        void loadFromFile (std::ifstream &file) override;
        void saveToFile (std::ofstream &file) const override;
    private:
        static const std::string M_NAME;
        static const GuiFabric<Cursor> m_cursorFabric;
        sf::Vertex m_points[6];
};

#endif