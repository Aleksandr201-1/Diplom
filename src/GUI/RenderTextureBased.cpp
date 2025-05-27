#include <GUI/RenderTextureBased.hpp>

//void RenderTextureBased::updateRender () {}

RenderTextureBased::RenderTextureBased () {
    m_clearColor = sf::Color::White;
    m_needRedraw = false;
    //setRenderSize(400, 400);
    //m_visibleContent.setTexture(m_render.getTexture());
}

RenderTextureBased::RenderTextureBased (uint64_t width, uint64_t height) : RenderTextureBased() {
    setRenderSize(width, height);
    //m_visibleContent.setTexture(m_render.getTexture(), true);
}

RenderTextureBased::RenderTextureBased (const sf::Vector2u &size) : RenderTextureBased() {
    setRenderSize(size);
    //m_visibleContent.setTexture(m_render.getTexture(), true);
}

RenderTextureBased::~RenderTextureBased () {}

void RenderTextureBased::setRenderSize (uint64_t width, uint64_t height) {
    setRenderSize(sf::Vector2u{width, height});
}

void RenderTextureBased::setRenderSize (const sf::Vector2u &size) {
    m_needRedraw = true;
    auto currSize = m_render.getSize();
    if (currSize.x < size.x || currSize.y < size.y) {
        currSize.x = std::max(size.x, currSize.x);
        currSize.y = std::max(size.y, currSize.y);
        m_render.create(currSize.x, currSize.y, true);
        m_visibleContent.setTexture(m_render.getTexture());
        m_visibleContent.setTextureRect({0, 0, size.x, size.y});
    }
}

void RenderTextureBased::setFonColor (const sf::Color &color) {
    m_needRedraw = true;
    m_clearColor = color;
}

void RenderTextureBased::loadFromFile (std::ifstream &file) {
    sf::Vector2u size;
    file.read(reinterpret_cast<char*>(&m_clearColor), sizeof(m_clearColor));
    file.read(reinterpret_cast<char*>(&m_needRedraw), sizeof(m_needRedraw));
    file.read(reinterpret_cast<char*>(&size), sizeof(size));
    if (size.x != 0 && size.y != 0) {
        m_render.create(size.x, size.y, true);
        m_visibleContent.setTexture(m_render.getTexture());
        m_visibleContent.setTextureRect({0, 0, size.x, size.y});
    }
}

void RenderTextureBased::saveToFile (std::ofstream &file) const {
    sf::Vector2u size = m_render.getSize();
    file.write(reinterpret_cast<const char*>(&m_clearColor), sizeof(m_clearColor));
    file.write(reinterpret_cast<const char*>(&m_needRedraw), sizeof(m_needRedraw));
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
}