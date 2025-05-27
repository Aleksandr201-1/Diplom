#ifndef INPUT_LINE_HPP
#define INPUT_LINE_HPP

#include <GUI/Interactable.hpp>
#include <GUI/TextBased.hpp>
#include <regex>

class InputLine : public Interactable, public TextBased {
    protected:
        void highLightRestart ();
        void replaceHighLightedWithString (const sf::String &string);
        void setCursorPos (uint64_t idx);
        void moveCursorToMouse (const sf::Vector2f &mousePos);
        void handleInput (const sf::Event &event);
        void keyPressedHandle (const sf::Event &event);
        void textEnteredHandle (const sf::Event &event);
        //void updateContent () override;
        void interact (const sf::Event &event, const sf::Vector2f &mousePos) override;
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
    public:
        InputLine ();
        InputLine (const sf::String &string);
        ~InputLine ();

        void setCheckRegex (const sf::String &regex);
        void setString (const sf::String &string) override;
        void copyToClipBoard () const;
        void setLimitOfSymbols (uint64_t limit);
        uint64_t getLimitOfSymbols () const;
        uint64_t getCursorPosition () const;
        sf::String getHighLightedText () const;
        double toFloat64 () const;
        int64_t toInt64 () const;
        float getCharPos (uint64_t idx);
        std::tuple<uint64_t, uint64_t> getCopyStartEnd ();

        using GUIElement::setSize;
        void setSize (const sf::Vector2f &size) override;
        void loadFromFile (std::ifstream &file) override;
        void saveToFile (std::ofstream &file) const override;
    private:
        static const std::string M_NAME;
        static const GuiFabric<InputLine> m_inputLineFabric;
        static const uint64_t M_LINE_OFFSET = 20;
        uint64_t m_copyStart, m_copyEnd;
        //sf::String m_checkRegex;
        //std::basic_regex<sf::Uint32> m_regex;
        std::wregex m_regex;
        //std::vector<sf::Uint32> m_forbiddenSymbols;
        sf::String m_highLightText;
        uint64_t m_limitOfSymbols, m_cursorPos;
};

#endif