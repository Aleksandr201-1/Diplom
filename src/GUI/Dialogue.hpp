#ifndef DIALOGUE_HPP
#define DIALOGUE_HPP

#include <GUI/ScrollingWindow.hpp>
#include <GUI/Interactable.hpp>
#include <GUI/TextBased.hpp>

// struct DialogueGeneral {
//     sf::Font font;
//     uint64_t size, speed;
//     sf::Color color;
    
//     DialogueGeneral (const std::string &pathToFont, uint64_t size = 12, uint64_t speed = 0, const sf::Color &color = sf::Color::Black);
//     ~DialogueGeneral ();
// };

class Dialogue : public Interactable, public TextBased {
    private:
        void updateContent () override;
        //взаимодействие
        void interact (const sf::Event &event, const sf::Vector2f &mousePos) override;
        //отрисовка
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
    public:
        Dialogue (const std::string &str);
        ~Dialogue ();
        void restart ();
        using GUIElement::setSize;
        void setSize (const sf::Vector2f &size) override;
        //void setTextSpeed (uint64_t s);
    private:
        //const std::string &str;
        sf::String m_string;
        bool isPlayed;
        sf::Clock clock;
        uint64_t m_currentPos, m_speed;
};

#endif