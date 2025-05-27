#ifndef DRAWABLE_HANDLER_HPP
#define DRAWABLE_HANDLER_HPP

#include <GUI/GUIElement.hpp>
#include <memory>

template <class T>
class DrawableHandler : public GUIElement {
    private:
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const;
    public:
        DrawableHandler ();
        DrawableHandler (std::unique_ptr<T> &&object);
        ~DrawableHandler ();

        void setObject (std::unique_ptr<T> &&object);
        std::unique_ptr<T> &getObject ();
        const std::unique_ptr<T> &getObject () const;
        using GUIElement::setSize;
        void setSize (const sf::Vector2f &size) override;
    private:
        std::unique_ptr<T> m_object;
};

template <class T>
void DrawableHandler<T>::drawVisible (sf::RenderTarget& target, sf::RenderStates states) const {
    states.transform *= getTransform();
    target.draw(*m_object.get(), states);
}

template <class T>
DrawableHandler<T>::DrawableHandler () {}

template <class T>
DrawableHandler<T>::DrawableHandler (std::unique_ptr<T> &&object) {
    setObject(object);
}

template <class T>
DrawableHandler<T>::~DrawableHandler () {}

template <class T>
void DrawableHandler<T>::setObject (std::unique_ptr<T> &&object) {
    m_object = std::move(object);
}

template <class T>
std::unique_ptr<T> &DrawableHandler<T>::getObject () {
    return m_object;
}

template <class T>
const std::unique_ptr<T> &DrawableHandler<T>::getObject () const {
    return m_object;
}

template <class T>
void DrawableHandler<T>::setSize (const sf::Vector2f &size) {
    m_size = size;
}

#endif