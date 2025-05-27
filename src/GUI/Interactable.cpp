#include <GUI/Interactable.hpp>
#include <iostream>

void Interactable::updateContent () {
    m_needUpdate = false;
}

void Interactable::clearSignals () {
    // for (uint64_t i = 0; i < m_signals.size(); ++i) {
    //     m_signals[i] = false;
    // }
    m_signals = 0;
}

// void Interactable::addSignal (Signal status) {
//     m_signals[static_cast<uint64_t>(status)] = true;
// }

void Interactable::updateFocusCondition (const sf::Event &event, const sf::Vector2f &mousePos) {
    //need trigger pressed
    //need trigger inside
    //need pressed inside to lose focus
    bool oldInteractFocus = m_interactFocus;
    if (m_trigger & to_underlying(ActivationTrigger::MOUSE)) {
        if (contains(mousePos)) { //inside
            handleSignal(Signal::TOUCHED);
            if (event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == m_mouse) { //mouse just pressed
                handleSignal(Signal::PRESSED);
                m_interactFocus = m_pressInsideToLoseFocus ? !m_interactFocus : true;
            } else if (sf::Mouse::isButtonPressed(m_mouse)) { //mouse alredy pressed
                //m_interactFocus = m_needTriggerInside && m_interactFocus ? true : m_interactFocus;
            } else if (event.type == sf::Event::MouseButtonReleased && event.mouseButton.button == m_mouse) { //mouse just released
                handleSignal(Signal::RELEASED);
                m_interactFocus = m_needTriggerPressed ? false : m_interactFocus;
            } else { //mouse already released
                //do nothing
                //m_interactFocus = m_needTriggerPressed ? m_interactFocus : false;
            }
        } else { //outside
            if (event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == m_mouse) { //mouse just pressed
                m_interactFocus = m_pressOutsideToLoseFocus ? false : m_interactFocus;
            } else if (sf::Mouse::isButtonPressed(m_mouse)) { //mouse alredy pressed
                //m_interactFocus = m_needTriggerPressed ? m_needTriggerInside ? false : m_interactFocus : m_interactFocus;
                //m_interactFocus = m_needTriggerInside && m_needTriggerPressed ? false : m_interactFocus;
            } else if (event.type == sf::Event::MouseButtonReleased && event.mouseButton.button == m_mouse) { //mouse just released
                m_interactFocus = m_needTriggerPressed ? false : m_interactFocus;
            } else { //mouse already released
                //do nothing
                //m_interactFocus = m_needTriggerPressed ? m_interactFocus : false;
            }
        }
    }
    if (m_trigger & to_underlying(ActivationTrigger::KEYBOARB)) {
        if (event.type == sf::Event::KeyPressed && event.key.code == m_key) {
            handleSignal(Signal::PRESSED);
            m_interactFocus = m_needTriggerPressed ? true : !m_interactFocus;
        } else {
            m_interactFocus = m_needTriggerPressed ? false : m_interactFocus;
        }
    }
    if (oldInteractFocus && !m_interactFocus) {
        handleSignal(Signal::LOST_FOCUS);
    }
    if (!oldInteractFocus && m_interactFocus) {
        handleSignal(Signal::RECEIVED_FOCUS);
    }
}

void Interactable::updateSizeCondition (const sf::Event &event, const sf::Vector2f &mousePos) {
    if (sf::Mouse::isButtonPressed(sf::Mouse::Left)) {

    }
}

Interactable::Interactable () {
    m_type = GUIType::INTERACTABLE;
    m_trigger = to_underlying(ActivationTrigger::MOUSE);
    m_mouse = sf::Mouse::Button::Left;
    m_key = sf::Keyboard::Enter;
    m_interactFocus = m_resizeFocus = false;
    m_pressInsideToLoseFocus = m_pressOutsideToLoseFocus = true;
    m_needTriggerPressed = m_needTriggerInside = false;
    m_needUpdate = true;
    m_signals = 0;
}

Interactable::~Interactable () {}

void Interactable::setInteractionFocus (bool focus) {
    m_interactFocus = focus;
}

bool Interactable::getInteractionFocus () const {
    return m_interactFocus;
}

void Interactable::setSignal (Signal status, const std::function<void (void)> &function) {
    m_workers[status] = function;
}

void Interactable::handleSignal (Signal signal) {
    //std::cout << "signal handling\n";
    uint32_t intVal = to_underlying(signal);
    size_t count = m_workers.count(signal);
    if (count == 1) {
        m_workers[signal]();
    } else if (count > 1) {
        std::cout << "Interactable::handleSignal: multiple signal handler detected";
        exit(1);
    }
    m_signals |= intVal;
}

void Interactable::needMousePressed (bool status) {
    m_needTriggerPressed = status;
}

bool Interactable::needUpdate () const {
    return m_needUpdate;
}

bool Interactable::getSignalStatus (Signal signal) const {
    return m_signals & to_underlying(signal);
}

// const std::array<bool, static_cast<uint64_t>(Interactable::Signal::COUNT_OF_SIGNALS)> &Interactable::getAllSignals () const {
//     return m_signals;
// }

uint32_t Interactable::getAllSignals () const {
    return m_signals;
}

void Interactable::setActivationTrigger (const sf::Mouse::Button &mouse) {
    m_mouse = mouse;
    m_trigger |= to_underlying(ActivationTrigger::MOUSE);
}

void Interactable::setActivationTrigger (const sf::Keyboard::Key &key) {
    m_key = key;
    m_trigger |= to_underlying(ActivationTrigger::KEYBOARB);
}

void Interactable::disableTrigger (ActivationTrigger trigger) {
    uint64_t intVal = to_underlying(trigger);
    if (m_trigger & intVal) {
        m_trigger -= intVal;
    }
}

void Interactable::handleEvent (const sf::Event &event, const sf::Vector2f &mousePos) {
    sf::Vector2f newMousePos = getTransform().getInverse().transformPoint(mousePos);
    updateFocusCondition(event, newMousePos);
    if (m_visible) {
        interact(event, newMousePos);
    }
    if (m_resizeFocus) {
        updateSizeCondition(event, newMousePos);
    }
}

void Interactable::update () {
    updateContent();
    clearSignals();
}

void Interactable::loadFromFile (std::ifstream &file) {
    file.read(reinterpret_cast<char*>(&m_trigger), sizeof(m_trigger));
    file.read(reinterpret_cast<char*>(&m_mouse), sizeof(m_mouse));
    file.read(reinterpret_cast<char*>(&m_key), sizeof(m_key));
    file.read(reinterpret_cast<char*>(&m_interactFocus), sizeof(m_interactFocus));
    file.read(reinterpret_cast<char*>(&m_resizeFocus), sizeof(m_resizeFocus));
    file.read(reinterpret_cast<char*>(&m_pressInsideToLoseFocus), sizeof(m_pressInsideToLoseFocus));
    file.read(reinterpret_cast<char*>(&m_pressOutsideToLoseFocus), sizeof(m_pressOutsideToLoseFocus));
    file.read(reinterpret_cast<char*>(&m_needTriggerPressed), sizeof(m_needTriggerPressed));
    file.read(reinterpret_cast<char*>(&m_needTriggerInside), sizeof(m_needTriggerInside));
    file.read(reinterpret_cast<char*>(&m_needUpdate), sizeof(m_needUpdate));
    file.read(reinterpret_cast<char*>(&m_signals), sizeof(m_signals));
    GUIElement::loadFromFile(file);
}

void Interactable::saveToFile (std::ofstream &file) const {
    file.write(reinterpret_cast<const char*>(&m_trigger), sizeof(m_trigger));
    file.write(reinterpret_cast<const char*>(&m_mouse), sizeof(m_mouse));
    file.write(reinterpret_cast<const char*>(&m_key), sizeof(m_key));
    file.write(reinterpret_cast<const char*>(&m_interactFocus), sizeof(m_interactFocus));
    file.write(reinterpret_cast<const char*>(&m_resizeFocus), sizeof(m_resizeFocus));
    file.write(reinterpret_cast<const char*>(&m_pressInsideToLoseFocus), sizeof(m_pressInsideToLoseFocus));
    file.write(reinterpret_cast<const char*>(&m_pressOutsideToLoseFocus), sizeof(m_pressOutsideToLoseFocus));
    file.write(reinterpret_cast<const char*>(&m_needTriggerPressed), sizeof(m_needTriggerPressed));
    file.write(reinterpret_cast<const char*>(&m_needTriggerInside), sizeof(m_needTriggerInside));
    file.write(reinterpret_cast<const char*>(&m_needUpdate), sizeof(m_needUpdate));
    file.write(reinterpret_cast<const char*>(&m_signals), sizeof(m_signals));
    GUIElement::saveToFile(file);
}