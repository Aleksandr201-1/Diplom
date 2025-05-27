#include <GUI/InputLine.hpp>

const std::string InputLine::M_NAME = "InputLine";
const InputLine::GuiFabric<InputLine> InputLine::m_inputLineFabric(InputLine::M_NAME);

void InputLine::highLightRestart () {
    m_copyStart = m_copyEnd = m_cursorPos;
}

void InputLine::replaceHighLightedWithString (const sf::String &string) {
    sf::String oldString = m_text.getString();
    uint64_t copyStart = std::min(m_copyStart, m_copyEnd);
    uint64_t copyEnd = std::max(m_copyStart, m_copyEnd);
    oldString.erase(copyStart, copyEnd - copyStart);
    oldString.insert(copyStart, string);
    m_cursorPos = string.getSize() + copyStart;
    highLightRestart();
    setCursorPos(m_cursorPos);
    m_text.setString(oldString);
    m_highLightText.clear();
    m_copyEnd = m_copyStart;
    handleSignal(Signal::CHANGED_VALUE);
}

void InputLine::setCursorPos (uint64_t idx) {
    m_cursorPos = idx;
    //m_cursor.setPosition({getCharPos(idx), 0});
}

void InputLine::moveCursorToMouse (const sf::Vector2f &mousePos) {
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

void InputLine::handleInput (const sf::Event &event) {
    if (event.type == sf::Event::KeyPressed) {
        keyPressedHandle(event);
        m_needUpdate = true;
    }
    if (event.type == sf::Event::TextEntered) {
        textEnteredHandle(event);
        m_needUpdate = true;
    }
}

void InputLine::keyPressedHandle (const sf::Event &event) {
    sf::String str = m_text.getString();
    float move = 0;
    uint64_t charSize = m_text.getCharacterSize();
    switch (event.key.code) {
        case sf::Keyboard::BackSpace:
            if (!m_highLightText.isEmpty()) {
                replaceHighLightedWithString("");
            } else if (m_cursorPos > 0) {
                --m_cursorPos;
                //glyth = m_text.getFont()->getGlyph(str[m_cursorPos], 20, false);
                move += m_text.getFont()->getGlyph(str[m_cursorPos], charSize, false).advance;
                //move += m_text.getLetterSpacing();
                if (m_cursorPos + 1 < str.getSize()) {
                    move += m_text.getFont()->getKerning(str[m_cursorPos], str[m_cursorPos + 1], charSize);
                }
                //m_cursor.move({-move, 0});
                str.erase(m_cursorPos, 1);
                m_text.setString(str);
                handleSignal(Signal::CHANGED_VALUE);
            }
            break;
        case sf::Keyboard::Delete:
            if (!m_highLightText.isEmpty()) {
                replaceHighLightedWithString("");
            } else if (m_cursorPos < str.getSize()) {
                str.erase(m_cursorPos, 1);
                m_text.setString(str);
                handleSignal(Signal::CHANGED_VALUE);
            }
            break;
        case sf::Keyboard::Enter:
            m_interactFocus = false;
            break;
        case sf::Keyboard::Right:
            if (m_copyStart != m_copyEnd) {
                m_cursorPos = std::max(m_copyStart, m_copyEnd);
                m_copyStart = m_copyEnd = 0;
                m_highLightText.clear();
            } else if (m_cursorPos < str.getSize()) {
                ++m_cursorPos;
                move += m_text.getFont()->getGlyph(str[m_cursorPos - 1], charSize, false).advance;
                //move += m_text.getLetterSpacing();
                if (m_cursorPos < str.getSize()) {
                    move += m_text.getFont()->getKerning(str[m_cursorPos - 1], str[m_cursorPos], charSize);
                }
                //m_cursor.move({move, 0});
            }
            break;
        case sf::Keyboard::Left:
            if (m_copyStart != m_copyEnd) {
                m_cursorPos = std::min(m_copyStart, m_copyEnd);
                m_copyStart = m_copyEnd = 0;
                m_highLightText.clear();
            } else if (m_cursorPos > 0) {
                --m_cursorPos;
                move += m_text.getFont()->getGlyph(str[m_cursorPos], charSize, false).advance;
                //move += m_text.getLetterSpacing();
                if (m_cursorPos + 1 < str.getSize()) {
                    move += m_text.getFont()->getKerning(str[m_cursorPos], str[m_cursorPos + 1], charSize);
                }
                //m_cursor.move({-move, 0});
            }
            break;
        default:
            break;
    }
}

void InputLine::textEnteredHandle (const sf::Event &event) {
    sf::String str = m_text.getString();
    float move = 0.f;
    uint64_t charSize = m_text.getCharacterSize();
    sf::String string;
    switch (event.text.unicode) {
        case 3: //Ctrl + C
            sf::Clipboard::setString(m_highLightText);
            break;
        case 22: //Ctrl + V
            string = sf::Clipboard::getString();
            if (!string.isEmpty()) {
                replaceHighLightedWithString(string);
            }
            handleSignal(Signal::CHANGED_VALUE);
            break;
        case '\t':
            break;
        case '\b':
            break;
        default:
            //if (std::regex_match(std::basic_string<wchar_t>(1, event.text.unicode), m_regex)) {
                if (!m_highLightText.isEmpty()) {
                    replaceHighLightedWithString(sf::String(event.text.unicode));
                } else if (str.getSize() + 1 < m_limitOfSymbols) {
                    str.insert(m_cursorPos, event.text.unicode);
                    ++m_cursorPos;
                    move += m_text.getFont()->getGlyph(str[m_cursorPos - 1], charSize, false).advance;
                    if (m_cursorPos < str.getSize()) {
                        move += m_text.getFont()->getKerning(str[m_cursorPos - 1], str[m_cursorPos], charSize);
                    }
                    //m_cursor.move({move, 0});
                    m_text.setString(str);
                }
                handleSignal(Signal::CHANGED_VALUE);
            //}
            break;
    }
}

void InputLine::interact (const sf::Event &event, const sf::Vector2f &mousePos) {
    if (m_interactFocus) {
        //mouse
        if (event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == m_mouse) {
            moveCursorToMouse(mousePos);
            m_copyStart = m_copyEnd = m_cursorPos;
            m_needUpdate = true;
        } else if (sf::Mouse::isButtonPressed(m_mouse) && event.type != sf::Event::LostFocus) {
            moveCursorToMouse(mousePos);
            m_copyEnd = m_cursorPos;
            //cursorRestart();
            if (m_copyStart != m_copyEnd) {
                uint64_t copyStart = std::min(m_copyStart, m_copyEnd);
                uint64_t copyEnd = std::max(m_copyStart, m_copyEnd);
                float startPos = getCharPos(copyStart), endPos = getCharPos(copyEnd);
                m_highLightText = m_text.getString().substring(copyStart, copyEnd - copyStart);
                handleSignal(Signal::SELECTED_TEXT);
            }
            m_needUpdate = true;
        }
        handleInput(event);
        m_size.x = m_text.getLocalBounds().width + M_LINE_OFFSET;
    }
}

void InputLine::drawVisible (sf::RenderTarget& target, sf::RenderStates states) const {
    target.draw(m_text, states);
}

InputLine::InputLine () : InputLine("") {}

InputLine::InputLine (const sf::String &string) : TextBased(string) {
    m_copyStart = m_copyEnd = 0;
    m_needTriggerPressed = false;
    m_needTriggerInside = false;
    m_pressInsideToLoseFocus = false;
    m_limitOfSymbols = 255;
    m_cursorPos = 0;
    m_size = {M_LINE_OFFSET, M_FONT_SIZE};
    m_needUpdate = true;
    m_needBorders = false;
    //sf::String str(L"*");
    //m_regex.assign(str.toWideString());
}

InputLine::~InputLine () {}

void InputLine::setCheckRegex (const sf::String &regex) {
    //m_regex.assign(regex.toWideString());
    //m_forbiddenSymbols = forbiddenSymbols;
}

void InputLine::setString (const sf::String &string) {
    TextBased::setString(string);
    m_copyStart = m_copyEnd = 0;
    m_cursorPos = 0;
    m_needUpdate = true;
}

void InputLine::copyToClipBoard () const {
    sf::Clipboard::setString(m_text.getString());
}

void InputLine::setLimitOfSymbols (uint64_t limit) {
    m_limitOfSymbols = limit;
}

uint64_t InputLine::getLimitOfSymbols () const {
    return m_limitOfSymbols;
}

uint64_t InputLine::getCursorPosition () const {
    return m_cursorPos;
}

sf::String InputLine::getHighLightedText () const {
    return m_highLightText;
}

double InputLine::toFloat64 () const {
    return std::atof(m_text.getString().toAnsiString().c_str());
}

int64_t InputLine::toInt64 () const {
    return std::atoll(m_text.getString().toAnsiString().c_str());
}

float InputLine::getCharPos (uint64_t idx) {
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

std::tuple<uint64_t, uint64_t> InputLine::getCopyStartEnd () {
    return std::make_tuple(m_copyStart, m_copyEnd);
}

void InputLine::setSize (const sf::Vector2f &size) {
    //m_needUpdate = true;
    //m_text.setCharacterSize(size.y);
    m_size = size;
}

void InputLine::loadFromFile (std::ifstream &file) {
    uint64_t size;
    //file.read(reinterpret_cast<char*>(&size), sizeof(size));
    //m_forbiddenSymbols.resize(size);
    //file.read(reinterpret_cast<char*>(&m_forbiddenSymbols[0]), sizeof(sf::Uint32) * size);
    file.read(reinterpret_cast<char*>(&m_limitOfSymbols), sizeof(m_limitOfSymbols));
    Interactable::loadFromFile(file);
    TextBased::loadFromFile(file);
}

void InputLine::saveToFile (std::ofstream &file) const {
    uint64_t size = M_NAME.size();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    file.write(&M_NAME[0], sizeof(char) * size);
    //size = m_forbiddenSymbols.size();
    //file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    //file.write(reinterpret_cast<const char*>(&m_forbiddenSymbols[0]), sizeof(sf::Uint32) * size);
    file.write(reinterpret_cast<const char*>(&m_limitOfSymbols), sizeof(m_limitOfSymbols));
    Interactable::saveToFile(file);
    TextBased::saveToFile(file);
}