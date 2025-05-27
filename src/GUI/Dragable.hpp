#ifndef DRAG_AND_DROP
#define DRAG_AND_DROP

#include <GUI/Interactable.hpp>

class Dragable : public Interactable {
    protected:
        void interact (const sf::Event &event, const sf::Vector2f &mousePos);
    public:
        Dragable ();
        virtual ~Dragable ();
    protected:
        sf::Vector2f m_dragableZone;
};

#endif