#include <GUI/MessageWindow.hpp>
#include <iostream>

const std::string MessageWindow::M_NAME = "MessageWindow";
const MessageWindow::GuiFabric<MessageWindow> MessageWindow::m_messageWindowFabric(MessageWindow::M_NAME);

void MessageWindow::recalculateSize () {
    //uint64_t charSize = m_text.getCharacterSize();
    //m_line.setPosition(m_offset, (static_cast<uint64_t>(m_size.y) - charSize) / 2);
    setRenderSize(m_size.x, m_size.y);
    m_visibleContent.setTextureRect({0, 0, m_size.x, m_size.y});
    m_button.setPosition((m_size.x - m_button.getSize().x) / 2, m_size.y - m_button.getSize().y - 10);
    m_needRedraw = true;
}

void MessageWindow::updateRender () {
    uint64_t charSize = m_text.getCharacterSize();
    m_render.clear(m_clearColor);
    sf::RectangleShape rect;
    sf::CircleShape circle;
    //sf::Vertex lines[2];
    float outlineThickness = 2.f;
    float rectIntend = 10.f;

    rect.setFillColor(sf::Color::Transparent);
    rect.setOutlineColor(sf::Color::Black);
    rect.setOutlineThickness(-outlineThickness);
    rect.setPosition(0, 0);
    rect.setSize({m_size.x, charSize + rectIntend});
    circle.setOutlineColor(sf::Color::Black);
    circle.setOutlineThickness(-2.f);
    circle.setRadius(static_cast<float>((charSize + rectIntend) / 2 - 2));
    circle.setPosition(m_size.x - circle.getRadius() * 2, 2);
    m_text.setOutlineColor(sf::Color::Black);
    m_text.setColor(sf::Color::White);
    m_text.setOutlineThickness(2.f);
    switch (m_messageType) {
        case Type::ERROR_MESS:
            m_text.setString(L"X");
            circle.setFillColor(sf::Color::Red);
            break;
        case Type::WARNING_MESS:
            m_text.setString(L"?");
            circle.setFillColor(sf::Color::Yellow);
            break;
        case Type::INFO_MESS:
            m_text.setString(L"!");
            circle.setFillColor(sf::Color::Blue);
            break;
        default:
            break;
    }

    m_render.draw(rect);
    float tmp = m_text.getFont()->getGlyph(m_text.getString()[0], charSize, false).advance;
    m_text.setPosition(m_size.x - circle.getRadius() - tmp / 2, 0);
    m_render.draw(circle);
    m_render.draw(m_text);

    m_text.setColor(sf::Color::Black);
    m_text.setOutlineThickness(0.f);

    m_text.setString(m_name);
    m_text.setPosition((m_size.x - m_text.getLocalBounds().width) / 2, 0);
    m_render.draw(m_text);

    rect.setSize({m_size.x, m_size.y - charSize + outlineThickness - 10});
    rect.setPosition(0, charSize + 10 - outlineThickness);
    m_render.draw(rect);

    m_text.setCharacterSize(15);
    m_text.setString(m_message);
    float ratio = m_text.getLocalBounds().width / (m_size.x - 10);
    uint64_t nLineCount = 1;
    if (ratio > 1.f) {
        nLineCount = (uint64_t)std::roundf(ratio + 0.5);
        sf::String str = m_text.getString();
        for (uint64_t k = 0; k < nLineCount - 1; ++k) {
            str.insert(str.getSize() / nLineCount + k, L"\n");
        }
        m_text.setString(str);
    }
    m_text.setPosition((m_size.x - m_text.getLocalBounds().width) / 2, 10 + rect.getSize().y / 2);
    m_render.draw(m_text);
    
    m_render.display();
    m_text.setCharacterSize(charSize);
    std::cout << "BE MessageWindow\n";
    m_needRedraw = false;
}

void MessageWindow::updateContent () {
    m_button.update();
    if (m_needRedraw) {
        updateRender();
        m_needRedraw = false;
    }
}

void MessageWindow::interact (const sf::Event &event, const sf::Vector2f &mousePos) {
    m_button.handleEvent(event, mousePos);
}

void MessageWindow::drawVisible (sf::RenderTarget& target, sf::RenderStates states) const {
    target.draw(m_visibleContent, states);
    target.draw(m_button, states);
}

MessageWindow::MessageWindow () : m_button(L"ОК") {
    setType(Type::INFO_MESS);
    setMessageName(L"Сообщение");
    setMessage(L"Текст сообщения");
    //m_button.setString(L"ОК");
    m_needRedraw = true;
}

MessageWindow::MessageWindow (Type type, const sf::String &name, const sf::String &message) : m_button(L"ОК") {
    setType(type);
    setMessageName(name);
    setMessage(message);
    //m_button.setString(L"ОК");
    m_needRedraw = true;
}

MessageWindow::~MessageWindow () {}

void MessageWindow::setType (Type type) {
    m_needUpdate = m_needRedraw = true;
    m_messageType = type;
}

void MessageWindow::setMessageName (const sf::String &name) {
    m_needUpdate = m_needRedraw = true;
    m_name = name;
}

void MessageWindow::setMessage (const sf::String &message) {
    m_needUpdate = m_needRedraw = true;
    m_message = message;
}

void MessageWindow::setSize (const sf::Vector2f &size) {
    m_size = size;
    recalculateSize();
}

void MessageWindow::loadFromFile (std::ifstream &file) {
    uint64_t size;
    std::basic_string<sf::Uint32> str;
    file.read(reinterpret_cast<char*>(&m_messageType), sizeof(m_messageType));

    file.read(reinterpret_cast<char*>(&size), sizeof(size));
    str.resize(size);
    file.read(reinterpret_cast<char*>(&str[0]), sizeof(str[0]) * size);
    m_message = str;

    file.read(reinterpret_cast<char*>(&size), sizeof(size));
    str.resize(size);
    file.read(reinterpret_cast<char*>(&str[0]), sizeof(str[0]) * size);
    m_name = str;

    m_button.loadFromFile(file);

    Interactable::loadFromFile(file);
    TextBased::loadFromFile(file);
    RenderTextureBased::loadFromFile(file);
}

void MessageWindow::saveToFile (std::ofstream &file) const {
    uint64_t size = M_NAME.size();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    file.write(&M_NAME[0], sizeof(char) * size);
    file.write(reinterpret_cast<const char*>(&m_messageType), sizeof(m_messageType));

    size = m_message.getSize();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    file.write(reinterpret_cast<const char*>(m_message.getData()), sizeof(m_message[0]) * size);

    size = m_name.getSize();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    file.write(reinterpret_cast<const char*>(m_name.getData()), sizeof(m_name[0]) * size);

    m_button.saveToFile(file);

    Interactable::saveToFile(file);
    TextBased::saveToFile(file);
    RenderTextureBased::saveToFile(file);
}