#include <GUI/Slider.hpp>

const std::string Slider::M_NAME = "Slider";
const Slider::GuiFabric<Slider> Slider::m_sliderFabric(Slider::M_NAME);

void Slider::interact (const sf::Event &event, const sf::Vector2f &mousePos) {
    float oldValue = m_value;
    float right = m_box.getPosition().x + m_box.getSize().x - m_slider.getSize().x - 5, left = m_box.getPosition().x + 5;
    float newPos = 0;
    float distance = right - left;
    float oneStep = distance / m_division;
    if (m_interactFocus) {
        newPos = mousePos.x - m_slider.getSize().x / 2;
        m_needUpdate = true;
    } else if (contains(mousePos) && event.type == sf::Event::MouseWheelMoved) {
        newPos = distance * m_value + left;
        if (m_needRound) {
            newPos -= float(event.mouseWheel.delta) * oneStep;
        } else {
            newPos -= float(event.mouseWheel.delta) * m_sensativity;
        }
        m_needUpdate = true;
    }
    if (m_needUpdate) {
        newPos = std::min(newPos, right);
        newPos = std::max(newPos, left);
        m_value = (newPos - left) / distance;
        if (m_needRound) {
            m_value = std::roundf(m_value * m_division) / m_division;
        }
        m_slider.setPosition(distance * m_value + left, 5);
        handleSignal(Signal::CHANGED_VALUE);
        handleSignal(Signal::VISUAL_CHANGE);
        m_needUpdate = false;
    }
}

void Slider::drawVisible (sf::RenderTarget& target, sf::RenderStates states) const {
    target.draw(m_box, states);
    target.draw(m_slider, states);
}

// Slider::Slider (const SpriteManager &manager, uint64_t first, uint64_t second, uint64_t min, uint64_t max): min(min), max(max) {
//     currState = (min + max) / 2;
//     needRound = true;
//     m_needMousePressed = true;
//     box = manager.at(first);
//     slider = manager.at(second);
//     slider.setTextureRect({0, 20, 20, 20});
//     size = {box.getGlobalBounds().width, box.getGlobalBounds().height};
// }

Slider::Slider () {
    m_needTriggerPressed = true;
    m_needRound = false;
    m_needUpdate = false;
    m_box.setOutlineColor(sf::Color::Black);
    m_box.setOutlineThickness(-2.f);
    //m_box.setFillColor(sf::Color::Transparent);
    m_box.setFillColor(sf::Color::White);
    m_box.setSize({400, 40});
    m_size = {400, 40};
    m_slider.setFillColor(sf::Color::Black);
    m_slider.setSize({20, 30});
    m_slider.setPosition(5, 5);
    m_division = 1;
    m_value = 0.f;
    m_min = 0.f;
    m_max = 1.f;
    m_sensativity = 10.f;
}

Slider::Slider (float min, float max, uint64_t division) : Slider() {
    setMinMax(min, max);
    m_division = division;
    if (m_division != 1) {
        m_needRound = true;
    } 
}

Slider::~Slider () {}

void Slider::setValue (float value) {
    m_value = value;
    float left = m_box.getPosition().x + m_box.getSize().x - m_slider.getSize().x - 5, right = m_box.getPosition().x + 5;
    float distance = left - right;
    if (m_needRound) {
        m_value = std::roundf(m_value * m_division) / m_division;
    }
    m_slider.setPosition(distance * m_value + right, 5);
    handleSignal(Signal::CHANGED_VALUE);
    handleSignal(Signal::VISUAL_CHANGE);
}

float Slider::getValue () const {
    return m_value * (m_max - m_min);
}

void Slider::setMinMax (float min, float max) {
    m_min = min;
    m_max = max;
    handleSignal(Signal::CHANGED_VALUE);
    handleSignal(Signal::VISUAL_CHANGE);
}

void Slider::setDivision (uint64_t division) {
    m_division = division;
}

void Slider::setSliderSize (float size) {
    m_slider.setSize({size, m_slider.getSize().y});
    handleSignal(Signal::VISUAL_CHANGE);
}

void Slider::setSensativity (float sensativity) {
    m_sensativity = sensativity;
}

void Slider::roundValue (bool status) {
    m_needRound = status;
}

void Slider::setSize (const sf::Vector2f &size) {
    m_size = size;
    //float ration = m_box.getSize().x / m_slider.getSize().x;
    m_box.setSize(m_size);
    //m_slider.setSize({m_size.x / ration, m_size.y - 10});
    m_slider.setSize({m_slider.getSize().x, m_size.y - 10});
    handleSignal(Signal::VISUAL_CHANGE);
}

void Slider::loadFromFile (std::ifstream &file) {
    file.read(reinterpret_cast<char*>(&m_needRound), sizeof(m_needRound));
    file.read(reinterpret_cast<char*>(&m_division), sizeof(m_division));
    file.read(reinterpret_cast<char*>(&m_value), sizeof(m_value));
    file.read(reinterpret_cast<char*>(&m_min), sizeof(m_min));
    file.read(reinterpret_cast<char*>(&m_max), sizeof(m_max));
    Interactable::loadFromFile(file);
    m_box.setSize(m_size);
    m_slider.setSize({m_slider.getSize().x, m_size.y - 10});
}

void Slider::saveToFile (std::ofstream &file) const {
    uint64_t size = M_NAME.size();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    file.write(&M_NAME[0], sizeof(char) * size);
    file.write(reinterpret_cast<const char*>(&m_needRound), sizeof(m_needRound));
    file.write(reinterpret_cast<const char*>(&m_division), sizeof(m_division));
    file.write(reinterpret_cast<const char*>(&m_value), sizeof(m_value));
    file.write(reinterpret_cast<const char*>(&m_min), sizeof(m_min));
    file.write(reinterpret_cast<const char*>(&m_max), sizeof(m_max));
    Interactable::saveToFile(file);
}