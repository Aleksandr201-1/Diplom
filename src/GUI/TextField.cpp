#include <GUI/TextField.hpp>
#include <iostream>

const std::string TextField::M_NAME = "TextField";
const TextField::GuiFabric<TextField> TextField::m_textFieldFabric(TextField::M_NAME);
const float TextField::M_CURSOR_X_OFFSET = 2.5f;

void TextField::recalculateSize () {
    uint64_t charSize = m_line.getText().getCharacterSize();
    m_line.setPosition(m_offset, (m_size.y - static_cast<float>(charSize)) / 2);
    setRenderSize(m_line.getLimitOfSymbols() * charSize, m_size.y);
    m_visibleContent.setTexture(m_render.getTexture());
    m_visibleContent.setTextureRect({m_offset, 0, m_size.x, m_size.y});
    m_needRedraw = true;
}

void TextField::updateRender () {
    m_render.clear(m_clearColor);
    m_render.draw(m_line);
    m_render.display();
    m_needRedraw = false;
    std::cout << "BE TF\n";
}

void TextField::cursorRestart () {
    m_clock.restart();
    m_visibleCursor = true;
}

void TextField::updateHighLightedText () {
    uint64_t charSize = m_line.getText().getCharacterSize();
    uint64_t copyStart, copyEnd;
    std::tie(copyStart, copyEnd) = m_line.getCopyStartEnd();
    if (copyStart != copyEnd) {
        m_visibleHighLight = true;
        m_copyStart = copyStart;
        m_copyEnd = copyEnd;
        copyStart = std::min(m_copyStart, m_copyEnd);
        copyEnd = std::max(m_copyStart, m_copyEnd);
        float startPos = m_line.getCharPos(copyStart), endPos = m_line.getCharPos(copyEnd);
        //m_highLightText.setString(m_text.getString().substring(copyStart, copyEnd - copyStart));
        if (startPos - m_offset < 0) {
            startPos = m_offset;
        }
        if (startPos - m_offset > m_size.x) {
            m_visibleHighLight = false;
        }
        if (endPos - m_offset < 0) {
            m_visibleHighLight = false;
        }
        if (endPos - m_offset > m_size.x) {
            endPos = m_size.x + m_offset;
        }
        m_blueBox.setPosition(startPos - m_offset + M_CURSOR_X_OFFSET, 0);
        m_blueBox.setSize({endPos - startPos, m_size.y});
        
        //handleSignal(Signal::SELECTED_TEXT);
    } else if (copyStart == copyEnd) {
        m_copyStart = copyStart;
        m_copyEnd = copyEnd;
        m_visibleHighLight = false;
    }
}

void TextField::updateCursorPos () {
    float move = -m_offset;
    float width = m_visibleContent.getTextureRect().width;
    const sf::Text &text = m_line.getText();
    const sf::String &str = text.getString();
    uint64_t charSize = text.getCharacterSize(), newPos = m_line.getCursorPosition();
    for (uint64_t i = 0; i < newPos; ++i) {
        move += text.getFont()->getGlyph(str[i], charSize, false).advance;
        if (i + 1 < str.getSize()) {
            move += text.getFont()->getKerning(str[i], str[i + 1], charSize);
        }
    }
    m_cursorPos = newPos;
    if (move < 0.f) {
        m_offset += move;
        //m_blueBox.move(move, 0);
        move = 0.f;
    }
    if (move > width) {
        m_offset += move - width;
        //m_blueBox.move(move - width, 0);
        move = width;
    }
    move += M_CURSOR_X_OFFSET;
    m_cursor.setPosition({move, 0});
    m_visibleContent.setTextureRect({m_offset, 0, m_size.x, m_size.y});
}

void TextField::updateContent () {
    if (m_clock.getElapsedTime().asSeconds() > 0.5f && m_interactFocus) {
        m_visibleCursor = !m_visibleCursor;
        m_clock.restart();
        handleSignal(Signal::VISUAL_CHANGE);
    }
    if (m_needUpdate) {
        m_line.update();
        updateCursorPos();
        updateHighLightedText();
        cursorRestart();
        m_needUpdate = false;
        handleSignal(Signal::VISUAL_CHANGE);
    }
    if (m_needRedraw) {
        updateRender();
    }
}

void TextField::interact (const sf::Event &event, const sf::Vector2f &mousePos) {
    if (event.type == sf::Event::MouseButtonPressed || event.type == sf::Event::KeyPressed) {
        cursorRestart();
    }
    if (m_interactFocus) {
        m_line.setInteractionFocus(true);
    } else {
        m_line.setInteractionFocus(false);
        m_visibleCursor = false;
    }
    uint64_t charSize = m_line.getText().getCharacterSize();
    //m_visibleCursor = m_interactFocus;
    m_line.handleEvent(event, mousePos + sf::Vector2f(m_offset, 0));
    m_needUpdate = m_line.needUpdate();
    if (m_line.getSignalStatus(Signal::CHANGED_VALUE)) {
        handleSignal(Signal::CHANGED_VALUE);
        handleSignal(Signal::VISUAL_CHANGE);
        if (!m_checker(m_line.getText().getString())) {
            handleSignal(Signal::ERROR_RECEIVED);
        }
        m_needRedraw = true;
    }
    if (getSignalStatus(Signal::LOST_FOCUS)) {
        handleSignal(Signal::VISUAL_CHANGE);
    }
    //m_needRedraw = m_line.getSignalStatus(Signal::CHANGED_VALUE);
    //m_needRedraw = m_needUpdate = m_interactFocus; //m_line.needUpdate();
}

void TextField::drawVisible (sf::RenderTarget& target, sf::RenderStates states) const {
    //target.draw(m_fon, states);
    target.draw(m_visibleContent, states);
    if (m_visibleHighLight && m_interactFocus) {
        target.draw(m_blueBox, states);
    }
    if (m_visibleCursor  && m_interactFocus) {
        target.draw(m_cursor, states);
    }
}

//void TextField::drawBorders (sf::RenderTarget& target, sf::RenderStates states) const {}

bool TextField::isFloat (const sf::String &string) {
    uint64_t pointCount = 0;
    for (uint64_t i = 0; i < string.getSize(); ++i) {
        if (string[i] == '.' || string[i] == ',') {
            ++pointCount;
            continue;
        }
        if (string[i] < '0' || string[i] > '9') {
            return false;
        }
    }
    return pointCount < 2;
}

bool TextField::isInt (const sf::String &string) {
    for (uint64_t i = 0; i < string.getSize(); ++i) {
        if (string[i] < '0' || string[i] > '9') {
            return false;
        }
    }
    return true;
}

TextField::TextField () : TextField("") {}

TextField::TextField (const sf::String &string) : m_line(string) {
    m_offset = 0;
    m_copyStart = m_copyEnd = 0;
    uint64_t charSize = m_line.getText().getCharacterSize();
    setRenderSize({m_line.getLimitOfSymbols() * charSize, charSize});
    m_visibleContent.setTextureRect({0, 0, 200, charSize});
    //m_fon.setFillColor(sf::Color::White);
    //m_fon.setOutlineColor(sf::Color::Black);
    //m_fon.setOutlineThickness(-2.f);
    //m_fon.setSize({100.f, m_line.getText().getCharacterSize()});
    m_blueBox.setFillColor(sf::Color(0, 0, 255, 127));
    m_checker = [] (const sf::String &string) {
        return true;
    };
    m_cursor.setFillColor(sf::Color::Black);
    m_cursor.setSize({2, m_line.getText().getCharacterSize()});
    m_cursor.move(M_CURSOR_X_OFFSET, 0);
    m_visibleCursor = false;
    m_pressInsideToLoseFocus = m_limitInputBySize = false;
    m_clock.restart();
    m_cursorPos = 0;
    m_needUpdate = m_needRedraw = true;
    m_size = {200, charSize};
}

TextField::~TextField () {}

InputLine &TextField::getInputLine () {
    m_needUpdate = m_needRedraw = true;
    return m_line;
}

const InputLine &TextField::getInputLine () const {
    return m_line;
}

void TextField::setCheckFunc (const std::function<bool (const sf::String &)> &checker) {
    m_checker = checker;
}

bool TextField::correctInput () const {
    return m_checker(m_line.getString());
}

void TextField::limitInputBySize (bool status) {
    m_limitInputBySize = status;
}

void TextField::setSize (const sf::Vector2f &size) {
    m_needUpdate = true;
    m_size = size;
    //m_line.setSize(size);
    //m_fon.setSize(size);
    m_cursor.setSize({m_cursor.getSize().x, size.y});
    recalculateSize();
}

void TextField::loadFromFile (std::ifstream &file) {
    Interactable::loadFromFile(file);
    RenderTextureBased::loadFromFile(file);
}

void TextField::saveToFile (std::ofstream &file) const {
    uint64_t size = M_NAME.size();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    file.write(&M_NAME[0], sizeof(char) * size);
    Interactable::saveToFile(file);
    RenderTextureBased::saveToFile(file);
}

// void TextField::update () {
//     if (m_clock.getElapsedTime().asSeconds() > 0.5f && m_interactFocus) {
//         m_visibleCursor = !m_visibleCursor;
//         m_clock.restart();
//     }
//     // if (m_needUpdate) {
//     //     updateCursorPos();
//     //     updateHighLightedText();
//     //     cursorRestart();
//     // }
//     Interactable::update();
//     if (m_needRedraw) {
//         updateRender();
//     }
// }