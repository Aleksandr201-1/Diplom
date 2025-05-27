#include <GUI/Legend.hpp>
#include <iostream>

const std::string Legend::M_NAME = "Legend";
const Legend::GuiFabric<Legend> Legend::m_legendFabric(Legend::M_NAME);

Legend::Record::Record () {}

Legend::Record::Record (const Line &line, const sf::String &string) : line(line), string(string) {}

Legend::Record::~Record () {}

void Legend::recalculateSize () {
    uint64_t maxIdx = 0;
    uint64_t letterSize = m_text.getCharacterSize();
    if (m_records.size() > 0) {
        for (uint64_t i = 0; i < m_records.size(); ++i) {
            if (m_records[i].string.getSize() > m_records[maxIdx].string.getSize()) {
                maxIdx = i;
            }
        }
        m_text.setString(m_records[maxIdx].string);
        m_size.x = 40.f + m_text.getGlobalBounds().width;
        m_size.y = 10.f + letterSize * m_records.size();
    } else {
        m_size.x = 40.f;
        m_size.y = 10.f;
    }
    setRenderSize(m_size.x, m_size.y);
}

void Legend::updateRender () {
    m_render.clear(m_clearColor);
    for (uint64_t i = 0; i < m_records.size(); ++i) {
        float letterSize = m_text.getCharacterSize();
        Line line = m_records[i].line;

        m_text.setString(m_records[i].string);
        m_text.setPosition(sf::Vector2f(30.f, letterSize * i));
        m_render.draw(m_text);

        line.position[0] = sf::Vector2f(5.f,  letterSize * 3 / 4 + letterSize * i);
        line.position[1] = sf::Vector2f(15.f, letterSize * 3 / 4 + letterSize * i);
        m_render.draw(line);

        line.position[0] = sf::Vector2f(15.f, letterSize * 3 / 4 + letterSize * i);
        line.position[1] = sf::Vector2f(25.f, letterSize * 3 / 4 + letterSize * i);
        m_render.draw(line);
    }
    m_render.display();
    m_needRedraw = false;
}

void Legend::updateContent () {
    if (m_needRedraw) {
        recalculateSize();
        updateRender();
    }
    m_needUpdate = false;
}

void Legend::interact (const sf::Event &event, const sf::Vector2f &mousePos) {}

void Legend::drawVisible (sf::RenderTarget& target, sf::RenderStates states) const {
    target.draw(m_visibleContent, states);
}

Legend::Legend () {
    m_needUpdate = m_needRedraw = true;
    m_text.setCharacterSize(15);
}

Legend::~Legend () {}

void Legend::add (const Record &record) {
    m_records.push_back(record);
    m_needUpdate = m_needRedraw = true;
}

void Legend::clear () {
    m_records.clear();
}

Legend::Record &Legend::operator[] (uint64_t i) {
    return m_records[i];
}

void Legend::setSize (const sf::Vector2f &size) {}

void Legend::loadFromFile (std::ifstream &file) {}

void Legend::saveToFile (std::ofstream &file) const {}

LegendAvailable::LegendAvailable () {}

LegendAvailable::~LegendAvailable () {}