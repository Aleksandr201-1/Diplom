#include <GUI/Selector.hpp>
#include <iostream>

const std::string Selector::M_NAME = "Selector";
const Selector::GuiFabric<Selector> Selector::m_selectorFabric(Selector::M_NAME);
const uint64_t Selector::CONTENT_CAP = 50;
const uint64_t Selector::MAX_LENGTH = 1000;
const uint64_t Selector::MIN_LENGTH = 50;

void Selector::recalculateSize () {
    float maxWidth = 0;
    float charSize = static_cast<float>(m_text.getCharacterSize());
    for (auto el : vars) {
        m_text.setString(el);
        auto size = m_text.getLocalBounds();
        maxWidth = std::max(maxWidth, size.width);
    }
    maxWidth = std::max(maxWidth + m_borderOffset * 2, static_cast<float>(Selector::MIN_LENGTH));
    maxWidth = std::min(maxWidth + m_borderOffset * 2, static_cast<float>(Selector::MAX_LENGTH));
    maxWidth = std::max(maxWidth, m_size.x);
    setRenderSize(maxWidth, charSize * Selector::CONTENT_CAP);
    m_currentContent.setTextureRect({0, 0, maxWidth, m_text.getCharacterSize()});
    m_visibleContent.setTextureRect({0, 0, maxWidth, m_text.getCharacterSize() * visibleOffset});
    checkRect.setSize({maxWidth, charSize});
    //visibleRect.setSize({maxWidth, charSize * visibleOffset});
    //currentRect.setSize({maxWidth, charSize});
    m_size = {maxWidth, charSize};
    handleSignal(Signal::VISUAL_CHANGE);
}

void Selector::updateRender () {
    float charSize = static_cast<float>(m_text.getCharacterSize());
    m_render.clear(m_clearColor);
    for (uint64_t i = 0; i < vars.size(); ++i) {
        m_text.setString(vars[i]);
        m_text.setPosition(m_borderOffset, charSize * i - charSize / 5);
        m_render.draw(m_text);
    }
    std::cout << "BE Slider\n";
    m_render.display();
    m_needRedraw = false;
}

void Selector::updateContent () {
    if (m_needRedraw) {
        recalculateSize();
        updateRender();
    }
    m_needUpdate = false;
}

void Selector::interact (const sf::Event &event, const sf::Vector2f &mousePos) {
    uint64_t charSize = m_text.getCharacterSize();
    int64_t offset = 0;
    int64_t oldStart = m_visibleContent.getTextureRect().top;
    uint64_t newIdx = -1;
    if (m_visibleContent.getGlobalBounds().contains(mousePos)) {
        if (m_interactFocus) {
            if (event.type == sf::Event::MouseWheelMoved) {
                float delta = event.mouseWheel.delta;
                offset = -(int64_t) delta * 5;
                if (oldStart + offset < 0) {
                    offset = -oldStart;
                } else if (oldStart + offset > (vars.size() - visibleOffset) * charSize) {
                    offset = (vars.size() - visibleOffset) * charSize - oldStart;
                }
                m_visibleContent.setTextureRect({0, oldStart + offset, m_size.x, charSize * visibleOffset});
                handleSignal(Signal::VISUAL_CHANGE);
            }
            newIdx = (static_cast<int64_t>(mousePos.y) + oldStart) / charSize;
            float newPos = static_cast<float>(charSize * newIdx) - oldStart - static_cast<float>(charSize);
            if (newPos < 0) {
                checkRect.setSize({m_size.x, charSize + newPos});
                newPos = 0;
            } else if (newPos > 2 * charSize) {
                checkRect.setSize({m_size.x, charSize + 2 * charSize - newPos});
            } else {
                checkRect.setSize({m_size.x, charSize});
            }
            checkRect.setPosition(0, newPos + charSize);
            handleSignal(Signal::VISUAL_CHANGE);
        }
        if (newIdx != -1) {
            m_checkIdx = newIdx - 1;
            handleSignal(Signal::VISUAL_CHANGE);
        }
        if (getSignalStatus(Signal::LOST_FOCUS)) {
            m_currIdx = m_checkIdx;
            //m_needUpdate = true;
            m_currentContent.setTextureRect({0, m_currIdx * charSize, m_size.x, charSize});
            handleSignal(Signal::CHANGED_VALUE);
            handleSignal(Signal::VISUAL_CHANGE);
        }
    }
    if (getSignalStatus(Signal::RECEIVED_FOCUS) || getSignalStatus(Signal::LOST_FOCUS)) {
        handleSignal(Signal::VISUAL_CHANGE);
    }
}

void Selector::drawVisible (sf::RenderTarget& target, sf::RenderStates states) const {
    if (m_interactFocus) {
        uint64_t charSize = m_text.getCharacterSize();
        target.draw(m_visibleContent, states);
        //target.draw(visibleRect, states);
        if (m_currIdx != -1) {
            target.draw(checkRect, states);
        }
        m_line.position[0] = {0, m_size.y};
        m_line.position[1] = {0, m_size.y + charSize * visibleOffset};
        target.draw(m_line, states);
        m_line.position[0] = {m_size.x, m_size.y + charSize * visibleOffset};
        target.draw(m_line, states);
        m_line.position[1] = {m_size.x, m_size.y};
        target.draw(m_line, states);
    }
    target.draw(m_currentContent, states);
    //target.draw(currentRect, states);
}

Selector::Selector () {
    m_countContentSize = false;
    m_borderOffset = 5.f;
    visibleOffset = 3;
    m_currentContent.setTexture(m_render.getTexture());
    m_currentContent.setTextureRect({0, 0, MAX_LENGTH, M_FONT_SIZE});
    m_visibleContent.setTextureRect({0, 0, MAX_LENGTH, M_FONT_SIZE * visibleOffset});
    m_currIdx = m_checkIdx = 0;

    checkRect.setFillColor(sf::Color(127, 127, 127, 127));
    checkRect.setOutlineColor(sf::Color::Black);
    checkRect.setOutlineThickness(-1.f);
    checkRect.setSize({MAX_LENGTH, (float)M_FONT_SIZE});

    // visibleRect.setFillColor(sf::Color::Transparent);
    // visibleRect.setOutlineColor(sf::Color::Black);
    // visibleRect.setOutlineThickness(-1.f);
    // visibleRect.setSize({MAX_LENGTH, M_FONT_SIZE * visibleOffset});

    // currentRect.setFillColor(sf::Color::Transparent);
    // currentRect.setOutlineColor(sf::Color::Black);
    // currentRect.setOutlineThickness(-1.f);
    // currentRect.setSize({MAX_LENGTH, M_FONT_SIZE});

    m_visibleContent.setPosition(0, (float)M_FONT_SIZE);
    //visibleRect.setPosition(0, (float)M_FONT_SIZE);

    m_size = {MIN_LENGTH, (float)M_FONT_SIZE};
}

Selector::Selector (const std::vector<sf::String> &list) : Selector() {
    setVarList(list);
}

Selector::~Selector () {}

void Selector::addVar (const sf::String &str) {
    m_needUpdate = m_needRedraw = true;
    vars.push_back(str);
    handleSignal(Signal::CHANGED_VALUE);
    handleSignal(Signal::VISUAL_CHANGE);
}

void Selector::setVarList (const std::vector<sf::String> &list) {
    m_needUpdate = m_needRedraw = true;
    vars = list;
    handleSignal(Signal::CHANGED_VALUE);
    handleSignal(Signal::VISUAL_CHANGE);
}

const sf::String &Selector::getCurrVal () const {
    return vars[m_currIdx];
}

uint64_t Selector::getCurrIdx () const {
    return m_currIdx;
}

void Selector::addContentToSize (bool status) {
    m_countContentSize = status;
}

void Selector::setVisibleOffset (uint64_t offset) {
    m_needUpdate = m_needRedraw = true;
    visibleOffset = offset;
}

uint64_t Selector::getVisibleOffset () const {
    return visibleOffset;
}

// void Selector::update () {
//     Interactable::update();
//     if (m_needRedraw) {
//         updateRender();
//     }
// }

void Selector::setSize (const sf::Vector2f &size) {
    recalculateSize();
    m_needRedraw = true;
    m_text.setCharacterSize(size.y);
    checkRect.setSize(size);
    //visibleRect.setSize({size.x, size.y * visibleOffset});
    //currentRect.setSize(size);
    m_currentContent.setTextureRect({0, 0, size.x, size.y});
    m_visibleContent.setTextureRect({0, 0, size.x, size.y * visibleOffset});
    m_visibleContent.setPosition(0, size.y);
    //visibleRect.setPosition(0, size.y);
    m_size = size;
}

void Selector::loadFromFile (std::ifstream &file) {
    uint64_t size;
    std::basic_string<sf::Uint32> str;
    file.read(reinterpret_cast<char*>(&m_borderOffset), sizeof(m_borderOffset));
    file.read(reinterpret_cast<char*>(&visibleOffset), sizeof(visibleOffset));
    file.read(reinterpret_cast<char*>(&m_currIdx), sizeof(m_currIdx));
    file.read(reinterpret_cast<char*>(&m_checkIdx), sizeof(m_checkIdx));
    file.read(reinterpret_cast<char*>(&size), sizeof(size));
    vars.resize(size);
    for (uint64_t i = 0; i < vars.size(); ++i) {
        file.read(reinterpret_cast<char*>(&size), sizeof(size));
        if (size > 0) {
            str.resize(size);
            file.read(reinterpret_cast<char*>(&str[0]), sizeof(str[0]) * size);
            vars[i] = str;
        }
    }
    Interactable::loadFromFile(file);
    TextBased::loadFromFile(file);
    RenderTextureBased::loadFromFile(file);
    recalculateSize();
}

void Selector::saveToFile (std::ofstream &file) const {
    uint64_t size = M_NAME.size();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    file.write(&M_NAME[0], sizeof(char) * size);
    file.write(reinterpret_cast<const char*>(&m_borderOffset), sizeof(m_borderOffset));
    file.write(reinterpret_cast<const char*>(&visibleOffset), sizeof(visibleOffset));
    file.write(reinterpret_cast<const char*>(&m_currIdx), sizeof(m_currIdx));
    file.write(reinterpret_cast<const char*>(&m_checkIdx), sizeof(m_checkIdx));
    size = vars.size();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (uint64_t i = 0; i < vars.size(); ++i) {
        size = vars[i].getSize();
        file.write(reinterpret_cast<const char*>(&size), sizeof(size));
        if (size > 0) {
            file.write(reinterpret_cast<const char*>(vars[i].getData()), sizeof(vars[i][0]) * size);
        }
    }
    Interactable::saveToFile(file);
    TextBased::saveToFile(file);
    RenderTextureBased::saveToFile(file);
}