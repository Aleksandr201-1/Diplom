#include <GUI/Button.hpp>

const std::string Button::M_NAME = "Button";
const Button::GuiFabric<Button> Button::m_buttonFabric(Button::M_NAME);
const sf::Color Button::M_REGULAR_COLOR = sf::Color(230, 230, 230), Button::M_TOUCH_COLOR = sf::Color(255, 255, 255), Button::M_CLICK_COLOR = sf::Color(200, 200, 200);

bool Button::touched (const sf::Event &event, const sf::Vector2f &mousePos) const {
    return button.getGlobalBounds().contains(mousePos) && event.type != sf::Event::LostFocus;
}

bool Button::clicked (const sf::Event &event, const sf::Vector2f &mousePos) const {
    return m_interactFocus;
}

void Button::changeShapeFromStatus () {
    switch (status) {
        case ButtonStatus::REGULAR:
            button.setFillColor(M_REGULAR_COLOR);
            break;
        case ButtonStatus::TOUCHED:
            button.setFillColor(M_TOUCH_COLOR);
            break;
        case ButtonStatus::CLICKED:
            button.setFillColor(M_CLICK_COLOR);
            break;
        default:
            break;
    }
}

void Button::updateContent () {
    // if (m_needUpdate) {
    //     button.setSize(m_size);
    //     m_text.setPosition({(m_size.x - m_text.getLocalBounds().width) / 2, 0});
    //     m_needUpdate = false;
    // }
}

void Button::interact (const sf::Event &event, const sf::Vector2f &mousePos) {
    ButtonStatus old = status;

    if (touched(event, mousePos)) {
        status = ButtonStatus::TOUCHED;
        handleSignal(Signal::TOUCHED);
    } else {
        status = ButtonStatus::REGULAR;
    }
    if (clicked(event, mousePos)) {
        status = ButtonStatus::CLICKED;
        handleSignal(Signal::PRESSED);
    }
    if (status != old) {
        m_needUpdate = true;
        changeShapeFromStatus();
        handleSignal(Signal::VISUAL_CHANGE);
    }
}

void Button::drawVisible (sf::RenderTarget& target, sf::RenderStates states) const {
    target.draw(button, states);
    target.draw(m_text, states);
}

Button::Button () {
    status = ButtonStatus::REGULAR;
    uint64_t charSize = m_text.getCharacterSize();
    button.setFillColor(sf::Color(230, 230, 230));
    button.setOutlineColor(sf::Color::Black);
    button.setOutlineThickness(-1.f);
    m_needTriggerPressed = true;
    m_needTriggerInside = true;
    m_pressOutsideToLoseFocus = false;
    m_size = {0, 0};
    m_needUpdate = false;
}

Button::Button (const sf::String &str) : Button() {
    m_text.setString(str);
    m_needUpdate = true;
    button.setSize({m_text.getLocalBounds().width + 20, 20 + 20});
    m_size = button.getSize();
    //m_text.setPosition({10, button.getSize().y / 2 - 10});
    m_text.setPosition({(m_size.x - m_text.getLocalBounds().width) / 2, (m_size.y - m_text.getLocalBounds().height) / 2});
}

Button::~Button () {}

//true, если на кнопку навели мышкой
bool Button::touched () const {
    return getSignalStatus(Signal::TOUCHED);
}

//true, если кнопку кликнули или ...
bool Button::pressed () const {
    return getSignalStatus(Signal::PRESSED);
}

void Button::setSize (const sf::Vector2f &size) {
    m_needUpdate = true;
    button.setSize(size);
    m_text.setPosition({(size.x - m_text.getLocalBounds().width) / 2, 0});
    m_size = size;
    handleSignal(Signal::VISUAL_CHANGE);
}

void Button::loadFromFile (std::ifstream &file) {
    file.read(reinterpret_cast<char*>(&status), sizeof(status));
    Interactable::loadFromFile(file);
    TextBased::loadFromFile(file);
    button.setSize(m_size);
}

void Button::saveToFile (std::ofstream &file) const {
    uint64_t size = M_NAME.size();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    file.write(&M_NAME[0], sizeof(char) * size);
    file.write(reinterpret_cast<const char*>(&status), sizeof(status));
    Interactable::saveToFile(file);
    TextBased::saveToFile(file);
}