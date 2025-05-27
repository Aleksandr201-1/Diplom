#ifndef MATH_TEXT_HPP
#define MATH_TEXT_HPP

#include <GUI/GUIElement.hpp>
#include <GUI/TextBased.hpp>
#include <GUI/RenderTextureBased.hpp>

class MathText : public GUIElement, public TextBased, public RenderTextureBased {
    private:
        enum NodeType {
            VALUE,
            FUNCTION
        };
        struct Node {
            NodeType type;
            std::string value;
        };
        using MathNode = std::unique_ptr<Node>;
    private:
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
    public:
        MathText ();
        MathText (const sf::String &string);
        ~MathText ();

        void setString (const sf::String &string);
    private:
        sf::String m_str;
        MathNode m_root;
};

#endif