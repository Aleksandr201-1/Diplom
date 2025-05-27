#include <GUI/Console.hpp>

Console::Command::Command () {}

Console::Command::Command (const sf::String &name, const std::function<void (const std::vector<sf::String> &)> &func) {
    m_name = name;
    m_func = func;
}

Console::Command::~Command () {}

bool Console::Command::operator== (const sf::String &str) {
    return str == m_name;
}

void Console::Command::operator() () {
    m_func({});
}

void Console::Command::operator() (const std::vector<sf::String> &args) {
    m_func(args);
}

void Console::updateRender () {
    m_render.clear(m_clearColor);
    uint64_t charSize = m_currentLine.getText().getCharacterSize();
    for (uint64_t i = 0; i < m_history.size(); ++i) {
        m_text.setString(m_userName + "$: " + m_history[i]);
        m_text.setPosition(0, (charSize + 5) * i);
    }
    m_text.setString(m_userName + "$: ");
    m_text.setPosition(0, (charSize + 5) * m_history.size());
    //m_currentLine.setPosition();
    m_render.display();
}

void Console::updateContent () {
    if (m_needUpdate) {
        m_currentLine.update();
        m_needUpdate = false;
    }
}

void Console::interact (const sf::Event &event, const sf::Vector2f &mousePos) {
    if (m_interactFocus) {
        
    }
}

void Console::drawVisible (sf::RenderTarget& target, sf::RenderStates states) const {
    target.draw(m_visibleContent, states);
    target.draw(m_currentLine, states);
}

Console::Console () {}

Console::~Console () {}

void Console::saveHistoryToFile (const std::string &file) {}

void Console::loadHistoryFromFile (const std::string &file) {}

std::vector<sf::String> Console::getHistory () {
    return m_history;
}

bool Console::receivedCommand () const {
    return getSignalStatus(Signal::RECEIVED_COMMAND);
}

sf::String Console::getLastCommand () {
    return m_history.back();
}

void Console::setHistoryLimit (uint64_t limit) {
    m_history.resize(limit);
    m_historyLimit = limit;
}

void Console::setLimitOfSymbols (uint64_t limit) {
    m_limitOfSymbols = limit;
}

void Console::setUserName (const sf::String &name) {
    m_needRedraw = true;
    m_userName = name;
}

sf::String Console::getUserName () const {
    return m_userName;
}

void Console::addCommand (const sf::String &name, const std::function<void (const std::vector<sf::String> &)> &func) {
    m_commands.push_back(Command(name, func));
    //m_commands[name] = func;
}

// void Console::addCommand (const Command &command) {
//     m_commands.push_back(command);
// }

// void Console::setCommandList (const std::vector<Command> &list) {
//     m_commands = list;
// }

void Console::clear () {
    m_history.clear();
    m_lineNumber = 0;
    m_needRedraw = true;
}

// void Console::update () {

// }

void Console::setSize (const sf::Vector2f &size) {}

std::ifstream &operator>> (std::ifstream &file, Console &console) {}

std::ofstream &operator<< (std::ofstream &file, const Console &console) {}