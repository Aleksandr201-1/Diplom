#include <GUI/GUIElement.hpp>

//std::map<std::string, std::string> GUIElement::m_fabric;
//std::vector<std::function<std::shared_ptr<GUIElement>()>> GUIElement::m_fab;

//sf::VertexArray GUIElement::m_line(sf::LineStrip, 4);

Line GUIElement::m_line(Line::Type::LINE, sf::Color::Black, 2.f);
//GUIElement::BorderInitializer GUIElement::m_borderInitializer;

void GUIElement::drawBorders (sf::RenderTarget& target, sf::RenderStates states) const {
    // m_line[0].position = {0, 0};
    // m_line[1].position = {m_size.x, 0};
    // m_line[2].position = {m_size.x, m_size.y};
    // m_line[3].position = {0, m_size.y};
    // target.draw(m_line, states);
    m_line.position[0] = {0, 0};
    m_line.position[1] = {m_size.x, 0};
    target.draw(m_line, states);
    m_line.position[0] = {m_size.x, m_size.y};
    target.draw(m_line, states);
    m_line.position[1] = {0, m_size.y};
    target.draw(m_line, states);
    m_line.position[0] = {0, 0};
    target.draw(m_line, states);
}

void GUIElement::draw (sf::RenderTarget& target, sf::RenderStates states) const {
    if (m_visible) {
        states.transform *= getTransform();
        drawVisible(target, states);
        if (m_needBorders) {
            drawBorders(target, states);
        }
    }
}

void GUIElement::printFabric () {
    for (auto &el : getFabric()) {
        std::cout << el.first << '\n';
    }
}

GUIElement::GUIElement () {
    m_type = GUIType::REGULAR;
    m_size = {0, 0};
    m_visible = m_resizable = m_needBorders = true;
}

GUIElement::~GUIElement () {}

void GUIElement::setVisibility (bool status) {
    m_visible = status;
}

bool GUIElement::isVisible () const {
    return m_visible;
}

bool GUIElement::isResizeable () const {
    return m_resizable;
}

bool GUIElement::needBorders () const {
    return m_needBorders;
}

void GUIElement::setBorders (bool status) {
    m_needBorders = status;
}

GUIElement::GUIType GUIElement::getType () const {
    return m_type;
}

sf::Vector2f GUIElement::getSize () const {
    return m_size;
}

void GUIElement::setSize (float width, float height) {
    setSize({width, height});
}

bool GUIElement::contains (const sf::Vector2f &pos) {
    return pos.x > 0.f && pos.x < m_size.x && pos.y > 0.f && pos.y < m_size.y;
}

void GUIElement::loadFromFile (const std::string &filename) {
    std::ifstream file(filename, std::ios_base::binary);
    loadFromFile(file);
    file.close();
}

void GUIElement::loadFromFile (std::ifstream &file) {
    //sf::Transform transform;
    sf::Vector2f position;
    float rotation;
    sf::Vector2f scale;
    sf::Vector2f origin;
    file.read(reinterpret_cast<char*>(&position), sizeof(position));
    file.read(reinterpret_cast<char*>(&rotation), sizeof(rotation));
    file.read(reinterpret_cast<char*>(&scale), sizeof(scale));
    file.read(reinterpret_cast<char*>(&origin), sizeof(origin));
    setPosition(position);
    setRotation(rotation);
    setScale(scale);
    setOrigin(origin);
    //file.read(reinterpret_cast<char*>(&transform), sizeof(transform));
    file.read(reinterpret_cast<char*>(&m_type), sizeof(m_type));
    file.read(reinterpret_cast<char*>(&m_size.x), sizeof(m_size.x));
    file.read(reinterpret_cast<char*>(&m_size.y), sizeof(m_size.y));
    file.read(reinterpret_cast<char*>(&m_visible), sizeof(m_visible));
    file.read(reinterpret_cast<char*>(&m_resizable), sizeof(m_resizable));
}

void GUIElement::saveToFile (const std::string &filename) const {
    std::ofstream file(filename, std::ios_base::binary);
    saveToFile(file);
    file.close();
}

void GUIElement::saveToFile (std::ofstream &file) const {
    //const sf::Transform &transform = getTransform();
    const sf::Vector2f &position = getPosition();
    float rotation = getRotation();
    const sf::Vector2f &scale = getScale();
    const sf::Vector2f &origin = getOrigin();
    file.write(reinterpret_cast<const char*>(&position), sizeof(position));
    file.write(reinterpret_cast<const char*>(&rotation), sizeof(rotation));
    file.write(reinterpret_cast<const char*>(&scale), sizeof(scale));
    file.write(reinterpret_cast<const char*>(&origin), sizeof(origin));
    //file.write(reinterpret_cast<const char*>(&transform), sizeof(transform));
    file.write(reinterpret_cast<const char*>(&m_type), sizeof(m_type));
    file.write(reinterpret_cast<const char*>(&m_size), sizeof(m_size));
    file.write(reinterpret_cast<const char*>(&m_visible), sizeof(m_visible));
    file.write(reinterpret_cast<const char*>(&m_resizable), sizeof(m_resizable));
}

std::map<std::string, std::function<std::shared_ptr<GUIElement>()>> &GUIElement::getFabric () {
    static std::map<std::string, std::function<std::shared_ptr<GUIElement>()>> m_fabric;
    return m_fabric;
}