#include <GUI/Dialogue.hpp>

// DialogueGeneral::DialogueGeneral (const std::string &pathToFont, uint64_t size, uint64_t speed, const sf::Color& color) : size(size), speed(speed),  color(color) {
//     if (!font.loadFromFile(pathToFont)) {
//         throw std::runtime_error("can't load font \"" + pathToFont + "\"");
//     }
// }

// DialogueGeneral::~DialogueGeneral () {}

// Dialogue::Dialogue (const DialogueGeneral &inf, const std::string &str) : general(inf) {
//     text.setFont(general.font);
//     text.setFillColor(general.color);
//     text.setCharacterSize(general.size);
//     text.setString(str);
//     currentChar = 0;
// }

void Dialogue::updateContent () {
    if (!isPlayed) {
        static uint64_t spc = 1000 / m_speed; // seconds per character
        int64_t time = clock.getElapsedTime().asMilliseconds();
        if (time > spc) {
            clock.restart();
            ++m_currentPos;
        }
        if (m_currentPos == m_string.getSize() - 1) {
            isPlayed = true;
        }
        m_text.setString(m_string.substring(0, m_currentPos + 1));
    }
}

void Dialogue::interact (const sf::Event &event, const sf::Vector2f &mousePos) {}

void Dialogue::drawVisible (sf::RenderTarget& target, sf::RenderStates states) const {
    target.draw(m_text, states);
}

Dialogue::Dialogue (const std::string &string) {
    setString(string);
    m_speed = 10;
}

Dialogue::~Dialogue () {}

void Dialogue::restart () {
    clock.restart();
    isPlayed = false;
    m_currentPos = 0;
}

void Dialogue::setSize (const sf::Vector2f &size) {}