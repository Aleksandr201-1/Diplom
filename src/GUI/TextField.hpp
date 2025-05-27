#ifndef TEXT_FIELD_HPP
#define TEXT_FIELD_HPP

#include <GUI/Interactable.hpp>
#include <GUI/RenderTextureBased.hpp>
#include <GUI/InputLine.hpp>
#include <functional>

class TextField : public Interactable, public RenderTextureBased {
    friend InputLine;
    private:
        void recalculateSize () override;
        void updateRender () override;
        void cursorRestart ();
        void updateHighLightedText ();
        void updateCursorPos ();
        void updateContent () override;
        void interact (const sf::Event &event, const sf::Vector2f &mousePos) override;
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
        //void drawBorders (sf::RenderTarget& target, sf::RenderStates states) const override;
    public:
        static bool isFloat (const sf::String &string);
        static bool isInt (const sf::String &string);
    public:
        TextField ();
        TextField (const sf::String &string);
        ~TextField ();

        InputLine &getInputLine ();
        const InputLine &getInputLine () const;
        void setCheckFunc (const std::function<bool (const sf::String &)> &checker);
        bool correctInput () const;
        void limitInputBySize (bool status);
        using GUIElement::setSize;
        void setSize (const sf::Vector2f &size) override;
        void loadFromFile (std::ifstream &file) override;
        void saveToFile (std::ofstream &file) const override;

        //void update () override;
    private:
        static const std::string M_NAME;
        static const GuiFabric<TextField> m_textFieldFabric;
        static const float M_CURSOR_X_OFFSET;
        float m_offset;
        uint64_t m_copyStart, m_copyEnd;
        InputLine m_line;
        sf::RectangleShape m_blueBox;
        bool m_visibleHighLight, m_limitInputBySize;
        std::function<bool (const sf::String &)> m_checker;
        sf::Clock m_clock;
        sf::RectangleShape m_cursor;
        bool m_visibleCursor;
        uint64_t m_cursorPos;
};

#endif