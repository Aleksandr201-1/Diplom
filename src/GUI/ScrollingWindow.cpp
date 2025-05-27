#include <GUI/ScrollingWindow.hpp>
#include <iostream>

const std::string ScrollingWindow::M_NAME = "ScrollingWindow";
const ScrollingWindow::GuiFabric<ScrollingWindow> ScrollingWindow::m_scrollingWindowFabric(ScrollingWindow::M_NAME);
const uint64_t ScrollingWindow::M_SLIDER_FREE_SPACE = 20;

void ScrollingWindow::recalculateSize () {
    std::cout << "recalculate\n";
    float maxWidth = 0, maxHeight = 0;
    for (const auto &el : m_content) {
        auto pos = el.second->getPosition();
        auto size = el.second->getSize();
        maxWidth = std::max(maxWidth, size.x + pos.x);
        maxHeight = std::max(maxHeight, size.y + pos.y);
    }
    if (m_xAxis) {
        maxHeight += m_sliderX.getSize().y;
    }
    if (m_yAxis) {
        maxWidth += m_sliderY.getSize().y;
    }
    //m_visibleContent.getGlobalBounds().width
    float prevWidth = static_cast<float>(m_render.getSize().x), prevHeight = static_cast<float>(m_render.getSize().y);
    if (maxWidth > prevWidth || maxHeight > prevHeight) {
        maxWidth = std::max(maxWidth, prevWidth);
        maxHeight = std::max(maxHeight, prevHeight);
        setRenderSize(maxWidth, maxHeight);
        m_availableOffset = {maxWidth, maxHeight};
        //m_visibleContent.setTextureRect({0, 0, m_size.x, m_size.y});
        m_needRedraw = true;
    }
    handleSignal(Signal::VISUAL_CHANGE);
    m_needRedraw = true;
    //m_needRecalculateSize = false;
    // maxWidth = std::max(maxWidth, static_cast<float>(m_render.getSize().x));
    // maxHeight = std::max(maxHeight, static_cast<float>(m_render.getSize().y));
    // m_availableOffset = {maxWidth, maxHeight};
    // setRenderSize(maxWidth, maxHeight);
    // m_visibleContent.setTextureRect({0, 0, m_size.x, m_size.y});
    //m_render.create(static_cast<uint64_t>(maxWidth),static_cast<uint64_t>(maxHeight));
    //m_size = {maxWidth, maxHeight};
}

void ScrollingWindow::updateRender () {
    m_render.clear(m_clearColor);
    for (auto &el : m_content) {
        m_render.draw(*el.second);
    }
    m_render.display();
    std::cout << "BE SW\n";
    m_needRedraw = false;
}

void ScrollingWindow::updateContent () {
    GUIContainer::updateContent();
    //getSignalStatus(Signal::VISUAL_CHANGE) || 
    if (getSignalStatus(Signal::RESIZED)) {
        recalculateSize();
    }
    if (m_needRedraw || getSignalStatus(Signal::VISUAL_CHANGE)) {
        updateRender();
    }
    m_sliderX.update();
    m_sliderY.update();
}

void ScrollingWindow::interact (const sf::Event &event, const sf::Vector2f &mousePos) {
    GUIContainer::interact(event, mousePos + sf::Vector2f(m_visibleContent.getTextureRect().left, m_visibleContent.getTextureRect().top));// + sf::Vector2f(m_visibleContent.getTextureRect().left, m_visibleContent.getTextureRect().top);
    if (m_xAxis) {
        m_sliderX.handleEvent(event, mousePos);
    }
    if (m_yAxis) {
        m_sliderY.handleEvent(event, mousePos);
    }
    if (m_sliderX.getSignalStatus(Signal::CHANGED_VALUE) || m_sliderY.getSignalStatus(Signal::CHANGED_VALUE)) {
        uint64_t distanceX = m_availableOffset.x - m_visibleContent.getTextureRect().width;
        uint64_t distanceY = m_availableOffset.y - m_visibleContent.getTextureRect().height;
        float xSliderValue = m_sliderX.getValue(), ySliderValue = m_sliderY.getValue();
        m_visibleContent.setTextureRect({distanceX * xSliderValue, distanceY * ySliderValue, m_size.x, m_size.y});
    }
    if (getSignalStatus(Signal::VISUAL_CHANGE)) {
        std::cout << "s";
        //handleSignal(Signal::VISUAL_CHANGE);
        //m_needRecalculateSize = true;
    }
    // if (m_availableOffset.x > m_visibleContent.getTextureRect().width && m_availableOffset.y > m_visibleContent.getTextureRect().height) {
    //     uint64_t distanceX = m_availableOffset.x - m_visibleContent.getTextureRect().width;
    //     uint64_t distanceY = m_availableOffset.y - m_visibleContent.getTextureRect().height;
    //     float xSliderValue = m_sliderX.getValue(), ySliderValue = m_sliderY.getValue();
    //     int64_t oldStart = m_visibleContent.getTextureRect().top;
    //     int64_t oldEnd = m_visibleContent.getTextureRect().height + oldStart;
    //     sf::IntRect rect;
    //     if (contains(mousePos)) {
    //         if (event.type == sf::Event::MouseWheelMoved && !m_interactFocus) {
    //             float offset = (-event.mouseWheel.delta * m_sensativity) / distanceX;
    //             ySliderValue += offset;
    //             if (ySliderValue < 0) {
    //                 ySliderValue = 0;
    //             }
    //             if (ySliderValue > 1.f) {
    //                 ySliderValue = 1.f;
    //             }
    //             m_sliderY.setValue(ySliderValue);
    //         }
    //     }
    //     rect.width = m_size.x;
    //     rect.height = m_size.y;
    //     //rect.left = std::min()
    //     m_visibleContent.setTextureRect({distanceX * xSliderValue, distanceY * ySliderValue, m_size.x, m_size.y});
    // }
}

void ScrollingWindow::drawVisible (sf::RenderTarget& target, sf::RenderStates states) const {
    target.draw(m_visibleContent, states);
    if (m_xAxis) {
        target.draw(m_sliderX, states);
    }
    if (m_yAxis) {
        target.draw(m_sliderY, states);
    }
}

ScrollingWindow::ScrollingWindow () {
    m_needRedraw = true;
    m_size = {400, 400};
    setRenderSize(static_cast<uint64_t>(m_size.x) * 2, static_cast<uint64_t>(m_size.y) * 2);
    //m_visibleContent.setTextureRect({0, 0, static_cast<uint64_t>(m_size.x), static_cast<uint64_t>(m_size.y)});
    m_availableOffset = {m_size.x * 2, m_size.y * 2};

    m_xAxis = m_yAxis = true;

    m_sliderX.setSize({m_size.x - M_SLIDER_FREE_SPACE + 2, M_SLIDER_FREE_SPACE});
    m_sliderX.setPosition(0, m_size.y - m_sliderX.getSize().y);
    m_sliderX.setSliderSize(30);
    m_sliderX.setValue(0);

    m_sliderY.rotate(90);
    m_sliderY.setSize({m_size.x, M_SLIDER_FREE_SPACE});
    m_sliderY.setPosition(m_size.x - m_sliderY.getSize().y + M_SLIDER_FREE_SPACE, 0);
    m_sliderY.setSliderSize(30);
}

ScrollingWindow::~ScrollingWindow () {}

void ScrollingWindow::setSensativity (float sensativity) {
    m_sliderX.setSensativity(sensativity);
    m_sliderX.setSensativity(sensativity);
}

void ScrollingWindow::setWindowSize (uint64_t width, uint64_t height) {
    setWindowSize({width, height});
}

void ScrollingWindow::setWindowSize (const sf::Vector2u &size) {

}

void ScrollingWindow::needXAxis (bool status) {
    m_xAxis = status;
}

void ScrollingWindow::needYAxis (bool status) {
    m_yAxis = status;
}

void ScrollingWindow::setSize (const sf::Vector2f &size) {
    float prevWidth = static_cast<float>(m_render.getSize().x), prevHeight = static_cast<float>(m_render.getSize().y);
    if (size.x > prevWidth || size.y > prevHeight) {
        float maxWidth = std::max(size.x, prevWidth);
        float maxHeight = std::max(size.y, prevHeight);
        setRenderSize(maxWidth, maxHeight);
        m_availableOffset = {maxWidth, maxHeight};
        m_needRedraw = true;
    }
    m_size = size;
    m_visibleContent.setTextureRect({0, 0, m_size.x, m_size.y});
    //m_needUpdate = true;
    //m_visibleContent.setTextureRect({0, 0, size.x, size.y});
    m_sliderX.setSize(size.x - M_SLIDER_FREE_SPACE + 2, m_sliderX.getSize().y);
    m_sliderX.setPosition(0, size.y - m_sliderX.getSize().y);
    //m_sliderY.setSize(m_sliderY.getSize().x, size.y);
    m_sliderY.setSize({size.y, m_sliderY.getSize().y});
    m_sliderY.setPosition(size.x - m_sliderY.getSize().y + M_SLIDER_FREE_SPACE, 0);
}

void ScrollingWindow::loadFromFile (std::ifstream &file) {
    uint64_t size;
    std::string str;
    file.read(reinterpret_cast<char*>(&size), sizeof(size));
    str.resize(size);
    file.read(&str[0], sizeof(char) * size);
    m_sliderX.loadFromFile(file);
    file.read(reinterpret_cast<char*>(&size), sizeof(size));
    str.resize(size);
    file.read(&str[0], sizeof(char) * size);
    m_sliderY.loadFromFile(file);
    file.read(reinterpret_cast<char*>(&m_xAxis), sizeof(m_xAxis));
    file.read(reinterpret_cast<char*>(&m_yAxis), sizeof(m_yAxis));
    file.read(reinterpret_cast<char*>(&m_availableOffset), sizeof(m_availableOffset));
    GUIContainer::loadFromFile(file);
    RenderTextureBased::loadFromFile(file);
}

void ScrollingWindow::saveToFile (std::ofstream &file) const {
    uint64_t size = M_NAME.size();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    file.write(&M_NAME[0], sizeof(char) * size);
    m_sliderX.saveToFile(file);
    m_sliderY.saveToFile(file);
    file.write(reinterpret_cast<const char*>(&m_xAxis), sizeof(m_xAxis));
    file.write(reinterpret_cast<const char*>(&m_yAxis), sizeof(m_yAxis));
    file.write(reinterpret_cast<const char*>(&m_availableOffset), sizeof(m_availableOffset));
    GUIContainer::saveToFile(file);
    RenderTextureBased::saveToFile(file);
}