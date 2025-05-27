#include <GUI/Dragable.hpp>

void Dragable::interact (const sf::Event &event, const sf::Vector2f &mousePos) {
    if (m_dragableZone) {
        
    }
    Interactable::interact(event, mousePos);
}

Dragable::Dragable () {}

Dragable::~Dragable () {}