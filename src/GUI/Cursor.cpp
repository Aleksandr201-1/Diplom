#include <GUI/Cursor.hpp>

const std::string Cursor::M_NAME = "Cursor";
const Cursor::GuiFabric<Cursor> Cursor::m_cursorFabric(Cursor::M_NAME);

void Cursor::updateContent () {}

void Cursor::interact (const sf::Event &event, const sf::Vector2f &mousePos) {
    if (event.type == sf::Event::MouseMoved) {
        m_points[0].position = sf::Vector2f(20, 20)   + mousePos;
        m_points[1].position = sf::Vector2f(0, 0)   + mousePos;
        m_points[2].position = sf::Vector2f(25, 12) + mousePos;
        m_points[3].position = sf::Vector2f(25, 25) + mousePos;
        m_points[4].position = sf::Vector2f(12, 25) + mousePos;
        m_points[5].position = sf::Vector2f(0, 0)   + mousePos;
    }
}

void Cursor::drawVisible (sf::RenderTarget& target, sf::RenderStates states) const {
    target.draw(m_points, 6, sf::TriangleFan, states);
}

Cursor::Cursor () {
    m_points[0].position = {20, 20};
    m_points[1].position = {0, 0};
    m_points[2].position = {25, 12};
    m_points[3].position = {25, 25};
    m_points[4].position = {12, 25};
    m_points[5].position = {0, 0};

    m_points[0].color = sf::Color(127, 127, 127);
    m_points[1].color = sf::Color::White;
    m_points[2].color = sf::Color::White;
    m_points[3].color = sf::Color::Black;
    m_points[4].color = sf::Color::Black;
    m_points[5].color = sf::Color(127, 127, 127);

    // for (uint64_t i = 0; i < 4; ++i) {
    //     m_points[i].color = sf::Color::Black;
    // }
}

Cursor::~Cursor () {}

void Cursor::setSize (const sf::Vector2f &size) {}

void Cursor::loadFromFile (std::ifstream &file) {
    file.read(reinterpret_cast<char*>(m_points), sizeof(sf::Vertex) * 6);
    Interactable::loadFromFile(file);
}

void Cursor::saveToFile (std::ofstream &file) const {
    uint64_t size = M_NAME.size();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    file.write(&M_NAME[0], sizeof(char) * size);
    file.write(reinterpret_cast<const char*>(m_points), sizeof(sf::Vertex) * 6);
    Interactable::saveToFile(file);
}