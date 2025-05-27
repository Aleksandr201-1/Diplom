#ifndef TEXT_BASED_HPP
#define TEXT_BASED_HPP

#include <SFML/Graphics.hpp>
#include <fstream>

class TextBased {
    public:
        TextBased ();
        TextBased (const sf::String &string);
        virtual ~TextBased ();

        void setCharacterSize (uint64_t charSize);
        void setFillColor (const sf::Color &color);
        void setFont (const sf::Font &font);
        const sf::Font &getUsedFont () const;
        virtual void setString (const sf::String &string);
        sf::String getString () const;
        sf::Text &getText ();
        const sf::Text &getText () const;

        virtual void loadFromFile (std::ifstream &file);
        virtual void saveToFile (std::ofstream &file) const;
    protected:
        static const uint64_t M_FONT_SIZE = 18;
        static const sf::Color M_FILL_COLOR;
        static sf::Font m_defaultFont;
        sf::Text m_text;
        //bool m_textChange;
    protected:
        static struct FontStaticConstructor {
            FontStaticConstructor () {
                m_defaultFont.loadFromFile("./source/consolas.ttf");
            }
        } m_fontStaticConsturctor;
};

#endif