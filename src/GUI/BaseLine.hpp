#ifndef BASE_LINE_HPP
#define BASE_LINE_HPP

#include <GUI/Interactable.hpp>
#include <GUI/TextBased.hpp>

class BaseLine : public Interactable, public TextBased {
    protected:
        void highLightRestart ();
        void setCursorPos (uint64_t idx);
        void moveCursorToMouse (const sf::Vector2f &mousePos);
        void updateContent () override;
        void interact (const sf::Event &event, const sf::Vector2f &mousePos) override;
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
    public:
        BaseLine ();
        BaseLine (const sf::String &string);
        virtual ~BaseLine ();

        sf::String find (const sf::String &pattern);
        void setString (const sf::String &string) override;
        void copyToClipBoard () const;
        uint64_t getCursorPosition () const;
        sf::String getHighLightedText () const;
        float getCharPos (uint64_t idx) const;
        std::tuple<uint64_t, uint64_t> getCopyStartEnd ();

        void loadFromFile (std::ifstream &file) override;
        void saveToFile (std::ofstream &file) const override;
    protected:
        static const uint64_t M_LINE_OFFSET = 20;
        uint64_t m_copyStart, m_copyEnd;
        //sf::String m_highLightText;
        uint64_t m_limitOfSymbols, m_cursorPos;
        sf::RectangleShape m_blueBox;
        bool m_visibleHighLight;
};

#endif