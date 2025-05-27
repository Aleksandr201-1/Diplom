#include <GUI/Bar.hpp>

const std::string Bar::M_NAME = "Bar";
const Bar::GuiFabric<Bar> Bar::M_BAR_FABRIC(Bar::M_NAME);

void Bar::updateContent () {
    float size = (m_param.get() > m_max) ? m_max : (m_param.get() < m_min) ? m_min : m_param.get();
    size /= m_max;
    if (m_rect.width * size != m_bar.getTextureRect().width) {
        m_bar.setTextureRect(sf::IntRect(m_rect.left, m_rect.top, m_rect.width * size, m_rect.height));
        handleSignal(Signal::VISUAL_CHANGE);
    }
}

void Bar::interact (const sf::Event &event, const sf::Vector2f &mousePos) {}

void Bar::drawVisible(sf::RenderTarget& target, sf::RenderStates states) const {
    target.draw(m_bar, states);
    target.draw(m_box, states);
}

Bar::Bar () : m_param(m_min), m_min(0), m_max(UINT64_MAX) {}

Bar::Bar (uint64_t &param, uint64_t min, uint64_t max) : m_param(param), m_min(min), m_max(max) {}

Bar::~Bar () {}

// void Bar::setSprites (const SpriteManager &manager, uint64_t first, uint64_t second) {
//     box = manager.at(first);
//     bar = manager.at(second);
//     m_rect = bar.getTextureRect();
// }

void Bar::setSprites (const sf::Sprite &box, const sf::Sprite &bar) {
    m_box = box;
    m_bar = bar;
    m_rect = m_bar.getTextureRect();
    handleSignal(Signal::VISUAL_CHANGE);
}

void Bar::setRef(uint64_t &param) {
    m_param = param;
    handleSignal(Signal::CHANGED_VALUE);
    handleSignal(Signal::VISUAL_CHANGE);
}

void Bar::setMinMax (uint64_t min, uint64_t max) {
    m_min = min;
    m_max = max;
}

void Bar::setSize (const sf::Vector2f &size) {
    scale(size.x / m_rect.width, size.y / m_rect.height);
    m_size = size;
    handleSignal(Signal::VISUAL_CHANGE);
}

void Bar::loadFromFile (std::ifstream &file) {
    uint64_t size;
    //file.read(reinterpret_cast<char*>(&size), sizeof(size));
    //file.read(&M_NAME[0], sizeof(char) * size);
    file.read(reinterpret_cast<char*>(&m_min), sizeof(m_min));
    file.read(reinterpret_cast<char*>(&m_max), sizeof(m_max));
    file.read(reinterpret_cast<char*>(&m_rect.left), sizeof(m_rect.left));
    file.read(reinterpret_cast<char*>(&m_rect.top), sizeof(m_rect.top));
    file.read(reinterpret_cast<char*>(&m_rect.width), sizeof(m_rect.width));
    file.read(reinterpret_cast<char*>(&m_rect.height), sizeof(m_rect.height));
    Interactable::loadFromFile(file);
}

void Bar::saveToFile (std::ofstream &file) const {
    uint64_t size = M_NAME.size();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    file.write(&M_NAME[0], sizeof(char) * size);
    file.write(reinterpret_cast<const char*>(&m_min), sizeof(m_min));
    file.write(reinterpret_cast<const char*>(&m_max), sizeof(m_max));
    file.write(reinterpret_cast<const char*>(&m_rect), sizeof(m_rect));
    Interactable::saveToFile(file);
    // file.write(reinterpret_cast<const char*>(&m_rect.left), sizeof(m_rect.left));
    // file.write(reinterpret_cast<const char*>(&m_rect.top), sizeof(m_rect.top));
    // file.write(reinterpret_cast<const char*>(&m_rect.width), sizeof(m_rect.width));
    // file.write(reinterpret_cast<const char*>(&m_rect.height), sizeof(m_rect.height));
}