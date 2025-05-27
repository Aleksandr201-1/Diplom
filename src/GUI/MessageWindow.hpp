#ifndef MESSAGE_WINDOW_HPP
#define MESSAGE_WINDOW_HPP

#include <GUI/Interactable.hpp>
#include <GUI/TextBased.hpp>
#include <GUI/RenderTextureBased.hpp>
#include <GUI/Button.hpp>
#include <cmath>

class MessageWindow : public Interactable, public TextBased, public RenderTextureBased {
    public:
        enum class Type : uint8_t {
            ERROR_MESS = 0,
            WARNING_MESS,
            INFO_MESS
        };
    private:
        void recalculateSize () override;
        void updateRender () override;
        void updateContent () override;
        void interact (const sf::Event &event, const sf::Vector2f &mousePos) override;
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
    public:
        MessageWindow ();
        MessageWindow (Type type, const sf::String &name, const sf::String &message);
        ~MessageWindow ();

        void setType (Type type);
        void setMessageName (const sf::String &name);
        void setMessage (const sf::String &message);

        using GUIElement::setSize;
        void setSize (const sf::Vector2f &size) override;
        void loadFromFile (std::ifstream &file) override;
        void saveToFile (std::ofstream &file) const override;
    private:
        static const std::string M_NAME;
        static const GuiFabric<MessageWindow> m_messageWindowFabric;
        Type m_messageType;
        //sf::RectangleShape m_rect;
        sf::String m_message, m_name;
        Button m_button;
};

#endif