#include <GUI/HeatMap.hpp>
#include <iostream>

void HeatMap::recalculateSize () {}

void HeatMap::updateRender () {
    // m_render.pushGLStates();
    // m_render.setActive(true);
    // glViewport(0, 0, m_render.getSize().x, m_render.getSize().y);
    // auto size = m_map.size();
    // float len = m_max - m_min;
    // if (len == 0.f) {
    //     len = 1.f;
    // }
    // glBegin(GL_QUADS);
    //     for (uint64_t i = 0; i < size.n - 1; ++i) {
    //         float curr;
    //         for (uint64_t j = 0; j < size.m - 1; ++j) {
    //             curr = (m_map(i, j) - m_min) / len;
    //             glColor3f(curr, (1.f - curr), 0.f);
    //             glVertex2f(j * m_squareSizeX, i * m_squareSizeY);
    //             curr = (m_map(i, j + 1) - m_min) / len;
    //             glColor3f(curr, (1.f - curr), 0.f);
    //             glVertex2f((j + 1) * m_squareSizeX, i * m_squareSizeY);
    //             curr = (m_map(i + 1, j + 1) - m_min) / len;
    //             glColor3f(curr, (1.f - curr), 0.f);
    //             glVertex2f((j + 1) * m_squareSizeX, (i + 1) * m_squareSizeY);
    //             curr = (m_map(i + 1, j) - m_min) / len;
    //             glColor3f(curr, (1.f - curr), 0.f);
    //             glVertex2f(j * m_squareSizeX, (i + 1) * m_squareSizeY);
    //         }
    //     }
    // glEnd();
    // m_render.setActive(false);
    // m_render.popGLStates();
    // // m_visibleContent.setTexture(m_render.getTexture());
    // // m_visibleContent.setTextureRect({0, 0, m_render.getSize().x , m_render.getSize().y});
    // m_render.display();
    m_needRedraw = false;
}

void HeatMap::updateContent () {
    if (m_needRedraw) {
        updateRender();
    }
    m_needUpdate = false;
}

void HeatMap::interact (const sf::Event &event, const sf::Vector2f &mousePos) {}

void HeatMap::drawVisible (sf::RenderTarget& target, sf::RenderStates states) const {
    // target.setActive(true);
    // glViewport(0, 0, target.getSize().x, target.getSize().y);
    // target.draw(m_visibleContent, states);
    target.pushGLStates();
    target.setActive(true);
    glViewport(0, 0, m_render.getSize().x, m_render.getSize().y);
    auto size = m_map.size();
    float len = m_max - m_min;
    if (len == 0.f) {
        len = 1.f;
    }
    glBegin(GL_QUADS);
        for (uint64_t i = 0; i < size.n - 1; ++i) {
            float curr;
            for (uint64_t j = 0; j < size.m - 1; ++j) {
                curr = (m_map(i, j) - m_min) / len;
                glColor3f(curr, (1.f - curr), 0.f);
                glVertex2f(j * m_squareSizeX, i * m_squareSizeY);
                curr = (m_map(i, j + 1) - m_min) / len;
                glColor3f(curr, (1.f - curr), 0.f);
                glVertex2f((j + 1) * m_squareSizeX, i * m_squareSizeY);
                curr = (m_map(i + 1, j + 1) - m_min) / len;
                glColor3f(curr, (1.f - curr), 0.f);
                glVertex2f((j + 1) * m_squareSizeX, (i + 1) * m_squareSizeY);
                curr = (m_map(i + 1, j) - m_min) / len;
                glColor3f(curr, (1.f - curr), 0.f);
                glVertex2f(j * m_squareSizeX, (i + 1) * m_squareSizeY);
            }
        }
    glEnd();
    target.popGLStates();
    target.setActive(false);
}

HeatMap::HeatMap () {
    m_min = INFINITY;
    m_max = -INFINITY;
    m_squareSizeX = m_squareSizeX = 400;
    m_size = {400, 400};
    setRenderSize(m_size.x, m_size.y);
}

HeatMap::~HeatMap () {}

void HeatMap::setMap (const Matrix<float> &map) {
    m_map = map;
    m_min = INFINITY;
    m_max = -INFINITY;
    for (auto el : m_map.toVector()) {
        m_min = std::min(m_min, el);
        m_max = std::max(m_max, el);
    }
    m_squareSizeX = m_size.x / m_map.size().m;
    m_squareSizeY = m_size.y / m_map.size().n;
}

void HeatMap::setSize (const sf::Vector2f &size) {}

void HeatMap::loadFromFile (std::ifstream &file) {}

void HeatMap::saveToFile (std::ofstream &file) const {}