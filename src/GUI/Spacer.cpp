#include <GUI/Spacer.hpp>

const std::string Spacer::M_NAME = "Spacer";
const Spacer::GuiFabric<Spacer> Spacer::m_spacerFabric(Spacer::M_NAME);

void Spacer::drawVisible (sf::RenderTarget& target, sf::RenderStates states) const {}

Spacer::Spacer () {
    setSize(0, 0);
}

Spacer::Spacer (float width, float height) {
    setSize(width, height);
}

Spacer::Spacer (const sf::Vector2f &size) {
    setSize(size);
}

Spacer::~Spacer () {}

void Spacer::setSize (const sf::Vector2f &size) {
    m_size = size;
}

void Spacer::saveToFile (std::ofstream &file) const {
    uint64_t size = M_NAME.size();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    file.write(&M_NAME[0], sizeof(char) * size);
}