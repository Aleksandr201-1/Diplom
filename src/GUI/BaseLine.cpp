#include <GUI/BaseLine.hpp>

void BaseLine::highLightRestart () {
    m_copyStart = m_copyEnd = m_cursorPos;
}

void BaseLine::setCursorPos (uint64_t idx) {
    m_cursorPos = idx;
    m_needUpdate = true;
}

void BaseLine::moveCursorToMouse (const sf::Vector2f &mousePos) {
    const sf::String &string = m_text.getString();
    uint64_t charSize = m_text.getCharacterSize(), stringSize = string.getSize();
    uint64_t newPos = 0;
    float width = 0, widthWithoutLastGlyph = 0, half = 0;
    if (string.getSize() == 0) {
        return;
    }
    while (width < mousePos.x && newPos < stringSize) {
        width += m_text.getFont()->getGlyph(string[newPos], charSize, false).advance;
        if (newPos + 1 < string.getSize()) {
            width += m_text.getFont()->getKerning(string[newPos], string[newPos + 1], charSize);
        }
        ++newPos;
    }
    if (newPos > 0) {
        widthWithoutLastGlyph = width - m_text.getFont()->getGlyph(string[newPos - 1], charSize, false).advance;
        widthWithoutLastGlyph -= m_text.getFont()->getKerning(string[newPos - 1], string[newPos], charSize);
        half = (width - widthWithoutLastGlyph) / 2;
        if (mousePos.x < widthWithoutLastGlyph + half) {
            width = widthWithoutLastGlyph;
            --newPos;
        }
    }
    //m_cursor.setPosition({width, 0});
    m_cursorPos = newPos;
}

void BaseLine::updateContent () {
    uint64_t charSize = getText().getCharacterSize();
    uint64_t copyStart, copyEnd;
    if (m_copyStart != m_copyEnd) {
        m_visibleHighLight = true;
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

void BaseLine::interact (const sf::Event &event, const sf::Vector2f &mousePos) {
    if (m_interactFocus) {
        //mouse
        if (event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == m_mouse) {
            moveCursorToMouse(mousePos);
            m_copyStart = m_copyEnd = m_cursorPos;
            m_needUpdate = true;
        } else if (sf::Mouse::isButtonPressed(m_mouse) && event.type != sf::Event::LostFocus) {
            moveCursorToMouse(mousePos);
            m_copyEnd = m_cursorPos;
            m_needUpdate = true;
        }
    }
}

void BaseLine::drawVisible (sf::RenderTarget& target, sf::RenderStates states) const {
    target.draw(m_text, states);
}

BaseLine::BaseLine () : BaseLine("") {}

BaseLine::BaseLine (const sf::String &string) : TextBased(string) {
    m_copyStart = m_copyEnd = 0;
    m_needTriggerPressed = false;
    m_needTriggerInside = false;
    m_pressInsideToLoseFocus = false;
    m_limitOfSymbols = 255;
    m_cursorPos = 0;
    m_size = {M_LINE_OFFSET, M_FONT_SIZE};
    m_needUpdate = true;
    m_needBorders = false;
}

BaseLine::~BaseLine () {}

void BaseLine::setString (const sf::String &string) {
    TextBased::setString(string);
    m_copyStart = m_copyEnd = 0;
    m_cursorPos = 0;
    m_needUpdate = true;
}

void BaseLine::copyToClipBoard () const {
    sf::Clipboard::setString(m_text.getString());
}

uint64_t BaseLine::getCursorPosition () const {
    return m_cursorPos;
}

sf::String BaseLine::getHighLightedText () const {
    sf::String highLightedText;
    if (m_copyStart != m_copyEnd) {
        uint64_t copyStart = std::min(m_copyStart, m_copyEnd);
        uint64_t copyEnd = std::max(m_copyStart, m_copyEnd);
        float startPos = getCharPos(copyStart), endPos = getCharPos(copyEnd);
        highLightedText = m_text.getString().substring(copyStart, copyEnd - copyStart);
    }
    return highLightedText;
}

float BaseLine::getCharPos (uint64_t idx) const {
    const sf::String &string = m_text.getString();
    uint64_t charSize = m_text.getCharacterSize(), stringSize = string.getSize();
    float pos = 0;
    for (uint64_t i = 0; i < idx; ++i) {
        pos += m_text.getFont()->getGlyph(string[i], charSize, false).advance;
        if (i + 1 < string.getSize()) {
            pos += m_text.getFont()->getKerning(string[i], string[i + 1], charSize);
        }
    }
    return pos;
}

std::tuple<uint64_t, uint64_t> BaseLine::getCopyStartEnd () {
    return std::make_tuple(m_copyStart, m_copyEnd);
}

void BaseLine::loadFromFile (std::ifstream &file) {
    file.read(reinterpret_cast<char*>(&m_copyStart), sizeof(m_copyStart));
    file.read(reinterpret_cast<char*>(&m_copyEnd), sizeof(m_copyEnd));
    file.read(reinterpret_cast<char*>(&m_limitOfSymbols), sizeof(m_limitOfSymbols));
    file.read(reinterpret_cast<char*>(&m_cursorPos), sizeof(m_cursorPos));
    file.read(reinterpret_cast<char*>(&m_visibleHighLight), sizeof(m_visibleHighLight));
    Interactable::loadFromFile(file);
    TextBased::loadFromFile(file);
}

void BaseLine::saveToFile (std::ofstream &file) const {
    decltype(m_copyStart) copyPos = 0;
    decltype(m_cursorPos) cursorPos = 0;
    decltype(m_visibleHighLight) visible = false;
    file.write(reinterpret_cast<const char*>(&copyPos), sizeof(copyPos));
    file.write(reinterpret_cast<const char*>(&copyPos), sizeof(copyPos));
    file.write(reinterpret_cast<const char*>(&m_limitOfSymbols), sizeof(m_limitOfSymbols));
    file.write(reinterpret_cast<const char*>(&cursorPos), sizeof(cursorPos));
    file.write(reinterpret_cast<const char*>(&visible), sizeof(visible));
    Interactable::saveToFile(file);
    TextBased::saveToFile(file);
}