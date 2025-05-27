#include <GUI/Label.hpp>

const std::string Label::M_NAME = "Label";
const Label::GuiFabric<Label> Label::m_labelFabric(Label::M_NAME);

void Label::drawVisible (sf::RenderTarget& target, sf::RenderStates states) const {
    target.draw(m_text, states);
}

Label::Label () : TextBased() {
    m_size = {0, 0};
    m_needBorders = false;
}

Label::Label (const sf::String &string) : TextBased(string) {
    m_size = {m_text.getLocalBounds().width + 5, m_text.getLocalBounds().height + 5};
}

Label::~Label () {}

// sf::Uint32 &Label::operator[] (uint64_t i) {
//     return m_text.getString()[i];
// }

sf::Uint32 Label::operator[] (uint64_t i) const {
    return m_text.getString()[i];
}

void Label::setSize (const sf::Vector2f &size) {
    //float xScale = size.x / m_size.x, yScale = size.y / m_size.y;
    //scale(xScale, yScale);
    m_size = size;
}

void Label::loadFromFile (std::ifstream &file) {
    GUIElement::loadFromFile(file);
    TextBased::loadFromFile(file);
}

void Label::saveToFile (std::ofstream &file) const {
    uint64_t size = M_NAME.size();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    file.write(&M_NAME[0], sizeof(char) * size);
    GUIElement::saveToFile(file);
    TextBased::saveToFile(file);
}