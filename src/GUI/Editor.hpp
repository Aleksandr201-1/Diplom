#ifndef EDITOR_HPP
#define EDITOR_HPP

#include <GUI/Interactable.hpp>
#include <GUI/TextBased.hpp>
#include <GUI/RenderTextureBased.hpp>
#include <GUI/InputLine.hpp>
#include <GUI/RichText.hpp>

class Editor : public Interactable, public RenderTextureBased {
    private:
        void updateRender () override;
        void updateContent () override;
        void interact (const sf::Event &event, const sf::Vector2f &mousePos) override;
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
    public:
        Editor ();
        ~Editor ();

        void setColorScheme (const std::map<sf::String, sf::Color> &colorScheme);
        void addToColorScheme (const sf::String &str, const sf::Color &color);
        void setAutocompletionRules (const std::map<sf::String, sf::String> &completionRules);
        void addToAutocompletionRules (const sf::String &str1, const sf::String &str2);
        void saveToFile (const std::string &file);
        bool lineNumbering (bool status);
        using GUIElement::setSize;
        void setSize (const sf::Vector2f &size) override;

        friend std::ostream &operator<< (std::ostream &output, const Editor &editor);
        friend std::istream &operator>> (std::istream &input, Editor &editor);
    private:
        bool m_lineNumbering;
        uint64_t m_intendLevel, m_intendSize, m_currentLine, m_currentPos, m_afterLineIndent;
        InputLine m_line;
        RichText m_text;
        std::map<sf::String, sf::Color> m_colorScheme;
        std::map<sf::String, sf::String> m_completionRules;
        std::vector<sf::String> m_content;
};

#endif