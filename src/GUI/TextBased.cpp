#include <GUI/TextBased.hpp>

const sf::Color TextBased::M_FILL_COLOR = sf::Color::Black;
sf::Font TextBased::m_defaultFont;
TextBased::FontStaticConstructor TextBased::m_fontStaticConsturctor;

TextBased::TextBased () : TextBased("") {}

TextBased::TextBased (const sf::String &string) {
    m_text.setFont(m_defaultFont);
    m_text.setFillColor(M_FILL_COLOR);
    m_text.setCharacterSize(M_FONT_SIZE);
    m_text.move(2.5, -(static_cast<float>(M_FONT_SIZE) / 5.f));
    setString(string);
}

TextBased::~TextBased () {}

void TextBased::setCharacterSize (uint64_t charSize) {
    //m_textChange = true;
    m_text.setCharacterSize(charSize);
}

void TextBased::setFillColor (const sf::Color &color) {
    //m_textChange = true;
    m_text.setFillColor(color);
}

void TextBased::setFont (const sf::Font &font) {
    //m_textChange = true;
    m_text.setFont(font);
}

const sf::Font &TextBased::getUsedFont () const {
    return *m_text.getFont();
}

void TextBased::setString (const sf::String &string) {
    //m_textChange = true;
    m_text.setString(string);
}

sf::String TextBased::getString () const {
    return m_text.getString();
}

sf::Text &TextBased::getText () {
    return m_text;
}

const sf::Text &TextBased::getText () const {
    return m_text;
}

void TextBased::loadFromFile (std::ifstream &file) {
    //sf::String str;
    std::basic_string<sf::Uint32> str;
    uint64_t size;
    sf::Color color;
    file.read(reinterpret_cast<char*>(&size), sizeof(size));
    if (size > 0) {
        str.resize(size);
        file.read(reinterpret_cast<char*>(&str[0]), sizeof(str[0]) * size);
        m_text.setString(str);
    }
    file.read(reinterpret_cast<char*>(&color), sizeof(color));
    file.read(reinterpret_cast<char*>(&size), sizeof(size));
    m_text.setCharacterSize(size);
}

void TextBased::saveToFile (std::ofstream &file) const {
    const sf::String &str = m_text.getString();
    uint64_t size = str.getSize();
    sf::Color color = m_text.getFillColor();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    if (size > 0) {
        file.write(reinterpret_cast<const char*>(str.getData()), sizeof(str[0]) * size);
    }
    file.write(reinterpret_cast<const char*>(&color), sizeof(color));
    size = m_text.getCharacterSize();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    //file.write(reinterpret_cast<const char*>(&m_maxWidth), sizeof(m_maxWidth));
}