#ifndef RICH_TEXT_HPP
#define RICH_TEXT_HPP

#include <GUI/GUIElement.hpp>

class RichText : public GUIElement {
    public:
        struct Style {
            const sf::Font *font;
            sf::Color color;
            uint64_t fontSize;
            uint32_t fontStyle;

            Style ();
            Style (const sf::Font &font, sf::Color color = sf::Color::Black, uint64_t fontSize = 28, uint32_t fontStyle = sf::Text::Style::Regular);
            Style (const Style &style);
            Style (Style &&style);
            ~Style ();

            Style &operator= (const Style &s);
            Style &operator= (Style &&s);
            friend bool operator== (const Style &s1, const Style &s2);
            friend bool operator!= (const Style &s1, const Style &s2);
        };
        struct StyleInterval {
            Style style;
            uint64_t begin, end;

            StyleInterval ();
            StyleInterval (const Style &style, uint64_t begin, uint64_t end);
            StyleInterval (const StyleInterval &style);
            StyleInterval (StyleInterval &&style);
            ~StyleInterval ();

            StyleInterval &operator= (const StyleInterval &s);
            StyleInterval &operator= (StyleInterval &&s);
            friend bool operator== (const StyleInterval &s1, const StyleInterval &s2);
            friend bool operator!= (const StyleInterval &s1, const StyleInterval &s2);
        };
    private:
        void calcMaxHeight ();
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
    public:
        RichText ();
        RichText (const sf::String &string, const Style &style);
        ~RichText ();

        void setFont (const sf::Font &font, uint64_t begin, uint64_t end);
        void setColor (const sf::Color &color, uint64_t begin, uint64_t end);
        void setFontSize (uint64_t fontSize, uint64_t begin, uint64_t end);
        void setFontStyle (uint32_t fontStyle, uint64_t begin, uint64_t end);
        void setStyleInterval (Style style, uint64_t begin, uint64_t end);
        void setUnderlineThickness (uint64_t thickness);
        void setUnderlineOffset (float offset);

        using GUIElement::setSize;
        void setSize (const sf::Vector2f &size) override;
        void loadFromFile (std::ifstream &file) override;
        void saveToFile (std::ofstream &file) const override;
    private:
        static const std::string M_NAME;
        static const GuiFabric<RichText> m_richTextFabric;
        mutable sf::Text m_text;
        sf::String m_string;
        float m_maxHeight;
        float m_lineOffset;
        uint64_t m_lineThickness;
        std::vector<StyleInterval> m_intervals;
};

#endif