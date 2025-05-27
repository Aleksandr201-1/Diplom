#ifndef GUI_STYLE_HPP
#define GUI_STYLE_HPP

#include <SFML/Graphics.hpp>
#include <memory>
#include <map>

class DrawableAndTransformable : public Drawable, public Transformable {};

class GuiStyle {
    public:
        GuiStyle ();
        ~GuiStyle ();

        std::shared_ptr<sf::Drawable> &operator[] (const std::string &id);
        const std::shared_ptr<sf::Drawable> &operator[] (const std::string &id) const;
    private:
        std::map<std::string, std::shared_ptr<DrawableAndTransformable>> m_style;
};

#endif