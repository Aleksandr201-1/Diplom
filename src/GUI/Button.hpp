#ifndef BUTTON_HPP
#define BUTTON_HPP

#include <GUI/Interactable.hpp>
#include <GUI/TextBased.hpp>

/*текстура хранит в себе "таблицу кнопок" в формате col.row
    _____________
    |1.1|1.2|1.3|
    |2.1|2.2|2.3|
    ...
    |n.1|n.2|n.3|
    ¯¯¯¯¯¯¯¯¯¯¯¯¯
    где:
    x.1 — обычное состояние (обязательно)
    x.2 — состояние при наведении мыши на кнопку (опционально)
    x.3 — состояние при нажатии на кнопку (опционально)
*/

class Button : public Interactable, public TextBased {
    public:
        enum class ButtonType {
            B_TOUCHABLE =   1 << 0,
            B_CLICKABLE =   1 << 1,
            B_TOUCH_SOUND = 1 << 2,
            B_CLICK_SOUND = 1 << 3,
            B_USE_TEXT =    1 << 4
        };

        enum class ButtonStatus {
            REGULAR,
            TOUCHED,
            CLICKED
        };
    private:
        //true, если на кнопку навели мышкой
        bool touched (const sf::Event &event, const sf::Vector2f& mousePos) const;
        //true, если кнопку кликнули
        bool clicked (const sf::Event &event, const sf::Vector2f& mousePos) const;
        void changeShapeFromStatus ();

        void updateContent () override;
        //взаимодействие
        void interact (const sf::Event &event, const sf::Vector2f &mousePos) override;
        //отрисовка
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
    public:
        //конструктор
        Button ();
        Button (const sf::String &str);
        //деструктор
        ~Button ();
        //true, если на кнопку навели мышкой
        bool touched () const;
        //true, если кнопку кликнули или ...
        bool pressed () const;
        using GUIElement::setSize;
        void setSize (const sf::Vector2f &size) override;
        void loadFromFile (std::ifstream &file) override;
        void saveToFile (std::ofstream &file) const override;
    private:
        static const std::string M_NAME;
        //elegant solution to fabric accidently found?!?
        static const GuiFabric<Button> m_buttonFabric;
        static const sf::Color M_REGULAR_COLOR, M_TOUCH_COLOR, M_CLICK_COLOR;
        ButtonStatus status;
        sf::RectangleShape button;
};

#endif