#ifndef BUTTON_BOX_HPP
#define BUTTON_BOX_HPP

#include <GUI/Interactable.hpp>
//#include <GUI/TextBased.hpp>
#include <GUI/Label.hpp>
#include <GUI/RenderTextureBased.hpp>

class ButtonBox : public Interactable, public RenderTextureBased {
    private:
        void recalculateSize ();
        void updateRender () override;
        void updateContent () override;
        void interact (const sf::Event &event, const sf::Vector2f &mousePos) override;
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
    public:
        ButtonBox ();
        ButtonBox (const std::vector<sf::String> content);
        ~ButtonBox ();

        void setItemList (const std::vector<sf::String> &strings);
        void addItem (const sf::String &string);
        uint64_t getCurrent () const;
        void setOffset (float offset);
        //void update () override;
        using GUIElement::setSize;
        void setSize (const sf::Vector2f &size) override;
        void loadFromFile (std::ifstream &file) override;
        void saveToFile (std::ofstream &file) const override;
    private:
        static const std::string M_NAME;
        static const GuiFabric<ButtonBox> m_buttonBoxFabric;
        float m_offset;
        uint64_t m_currentId;
        std::vector<Label> m_strings;
        sf::CircleShape m_currentShape, m_normalShape;
};

#endif