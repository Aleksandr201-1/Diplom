#include <GUI/ButtonBox.hpp>
#include <iostream>

const std::string ButtonBox::M_NAME = "ButtonBox";
const ButtonBox::GuiFabric<ButtonBox> ButtonBox::m_buttonBoxFabric(ButtonBox::M_NAME);

void ButtonBox::recalculateSize () {
    if (m_strings.empty()) {
        return;
    }
    float maxWidth = 0, maxHeight = 0;
    for (uint64_t i = 0; i < m_strings.size(); ++i) {
        maxWidth = std::max(maxWidth, m_strings[i].getSize().x + m_normalShape.getRadius() * 2 + 5);
        maxHeight += m_strings[i].getSize().y + m_offset;
    }
    maxHeight -= m_offset;
    m_size = {maxWidth, maxHeight};
    setRenderSize(static_cast<uint64_t>(m_size.x), static_cast<uint64_t>(m_size.y));
}

void ButtonBox::updateRender () {
    if (m_strings.empty()) {
        return;
    }
    m_render.clear(m_clearColor);
    float yLabelSize = m_strings[0].getSize().y;
    float circleSize = m_normalShape.getRadius() * 2;
    for (uint64_t i = 0; i < m_strings.size(); ++i) {
        m_normalShape.setPosition(0, i * (m_offset + yLabelSize) + (yLabelSize - circleSize) / 2);
        m_strings[i].setPosition(circleSize + 5, i * (m_offset + yLabelSize));
        m_render.draw(m_normalShape);
        m_render.draw(m_strings[i]);
    }
    m_render.display();
    //m_visibleContent.setTexture(m_render.getTexture());
    m_needRedraw = false;
}

void ButtonBox::updateContent () {
    if (m_needUpdate) {
        recalculateSize();
        m_needUpdate = false;
    }
    if (m_needRedraw) {
        updateRender();
        m_needRedraw = false;
    }
    //recalculateSize();
    //setRenderSize(static_cast<uint64_t>(m_size.x), static_cast<uint64_t>(m_size.y));
}

void ButtonBox::interact (const sf::Event &event, const sf::Vector2f &mousePos) {
    float yLabelSize = m_strings[0].getSize().y;
    float yMousePos = mousePos.y;
    float circleSize = m_normalShape.getRadius() * 2;
    float currentCircleSize = m_currentShape.getRadius() * 2;
    if (contains(mousePos) && event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == m_mouse) {
        for (uint64_t i = 0; i < m_strings.size(); ++i) {
            float curr = i * (m_offset + yLabelSize);
            if (yMousePos > curr && yMousePos < curr + yLabelSize) {
                m_currentId = i;
                m_currentShape.setPosition((circleSize - currentCircleSize) / 2, i * (m_offset + yLabelSize) + (yLabelSize - currentCircleSize) / 2);
                handleSignal(Signal::CHANGED_VALUE);
                break;
            }
        }
    }
    m_interactFocus = false;
}

void ButtonBox::drawVisible (sf::RenderTarget& target, sf::RenderStates states) const {
    target.draw(m_visibleContent, states);
    if (m_currentId != -1) {
        target.draw(m_currentShape, states);
    }
}

ButtonBox::ButtonBox () {
    //setRenderSize(200, 400);
    //m_size = {200, 400};
    m_offset = 20;
    m_currentShape.setRadius(7.f);
    m_currentShape.setFillColor(sf::Color::Black);
    m_currentShape.setOutlineColor(sf::Color::Black);
    m_normalShape.setRadius(10.f);
    m_normalShape.setFillColor(sf::Color::White);
    m_normalShape.setOutlineColor(sf::Color::Black);
    m_normalShape.setOutlineThickness(-1.f);
    m_currentId = -1;
}

ButtonBox::ButtonBox (const std::vector<sf::String> content) : ButtonBox() {
    setItemList(content);
}

ButtonBox::~ButtonBox () {}

void ButtonBox::setItemList (const std::vector<sf::String> &strings) {
    m_needRedraw = m_needUpdate = true;
    m_currentId = 0;
    m_strings.clear();
    for (uint64_t i = 0; i < strings.size(); ++i) {
        m_strings.push_back(std::move(Label(strings[i])));
    }
    float yLabelSize = m_strings[0].getSize().y;
    float circleSize = m_normalShape.getRadius() * 2;
    float currentCircleSize = m_currentShape.getRadius() * 2;
    m_currentShape.setPosition((circleSize - currentCircleSize) / 2, (yLabelSize - currentCircleSize) / 2);
}

void ButtonBox::addItem (const sf::String &string) {
    m_needRedraw = m_needUpdate = true;
    m_strings.push_back(Label(string));
}

uint64_t ButtonBox::getCurrent () const {
    return m_currentId;
}

void ButtonBox::setOffset (float offset) {
    m_needUpdate = m_needRedraw = true;
    m_offset = offset;
}

// void ButtonBox::update () {
//     if (m_needUpdate) {
//         recalculateSize();
//     }
//     Interactable::update();
//     if (m_needRedraw) {
//         updateRender();
//     }
// }

void ButtonBox::setSize (const sf::Vector2f &size) {
    recalculateSize();
    m_needRedraw = true;
    uint64_t count = m_strings.size() - 1;
    float newOffset = (size.y - m_size.y + count * m_offset) / count;
    std::cout << "new offset: " << newOffset << "\n";
    if (newOffset > 0.f) {
        m_size = size;
        m_offset = newOffset;
    }
}

void ButtonBox::loadFromFile (std::ifstream &file) {
    uint64_t size;
    std::string str;
    file.read(reinterpret_cast<char*>(&m_offset), sizeof(m_offset));
    file.read(reinterpret_cast<char*>(&m_currentId), sizeof(m_currentId));
    file.read(reinterpret_cast<char*>(&size), sizeof(size));
    m_strings.resize(size);
    for (uint64_t i = 0; i < size; ++i) {
        file.read(reinterpret_cast<char*>(&size), sizeof(size));
        str.resize(size);
        file.read(&str[0], sizeof(char) * size);
        m_strings[i].loadFromFile(file);
    }
    Interactable::loadFromFile(file);
    RenderTextureBased::loadFromFile(file);
}

void ButtonBox::saveToFile (std::ofstream &file) const {
    uint64_t size = M_NAME.size();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    file.write(&M_NAME[0], sizeof(char) * size);
    file.write(reinterpret_cast<const char*>(&m_offset), sizeof(m_offset));
    file.write(reinterpret_cast<const char*>(&m_currentId), sizeof(m_currentId));
    size = m_strings.size();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (const auto &label : m_strings) {
        label.saveToFile(file);
        // size = label.getString().getSize();
        // file.write(reinterpret_cast<const char*>(&size), sizeof(size));
        // file.write(reinterpret_cast<const char*>(label.getData()), sizeof(label[0]) * size);
        // file.write(&el.getString()[0], sizeof(char) * size);
    }
    Interactable::saveToFile(file);
    RenderTextureBased::saveToFile(file);
}