#include <GUI/Page.hpp>

//const std::string Page::M_NAME = "Page";

void Page::recalculateSize () {
    float maxWidth = 0, maxHeight = 0;
    for (const auto &el : m_content) {
        auto pos = el.second->getPosition();
        auto size = el.second->getSize();
        maxWidth = std::max(maxWidth, size.x + pos.x);
        maxHeight = std::max(maxHeight, size.y + pos.y);
    }
    m_size = {maxWidth, maxHeight};
}

Page::Page () {}

Page::~Page () {}

void Page::setSize (const sf::Vector2f &size) {
    m_size = size;
}

// void Page::saveToFile (std::ofstream &file) const {
//     //uint64_t size = M_NAME.size();
//     //file.write(reinterpret_cast<const char*>(&size), sizeof(size));
//     //file.write(&M_NAME[0], sizeof(char) * size);
//     GUIContainer::saveToFile(file);
// }