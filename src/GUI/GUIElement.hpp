#ifndef GUI_ELEMENT_HPP
#define GUI_ELEMENT_HPP

#include <SFML/Graphics.hpp>
#include <SFExtensions/Line.hpp>
#include <fstream>
#include <memory>
#include <map>
#include <functional>
#include <iostream>

class GUIElement : public sf::Transformable, public sf::Drawable {
    public:
        enum class GUIType : uint64_t {
            REGULAR = 0,
            INTERACTABLE,
            CONTAINER
        };
    protected:
        virtual void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const = 0;
        virtual void drawBorders (sf::RenderTarget& target, sf::RenderStates states) const;
        void draw (sf::RenderTarget& target, sf::RenderStates states) const override;
    public:
        //to delete
        static void printFabric ();
        GUIElement ();
        virtual ~GUIElement ();
        void setVisibility (bool status);
        bool isVisible () const;
        bool isResizeable () const;
        bool needBorders () const;
        void setBorders (bool status);
        GUIType getType () const;
        void setSize (float width, float height);
        virtual void setSize (const sf::Vector2f &v) = 0;
        sf::Vector2f getSize () const;
        bool contains (const sf::Vector2f &pos);
        void loadFromFile (const std::string &filename);
        virtual void loadFromFile (std::ifstream &file);
        void saveToFile (const std::string &filename) const;
        virtual void saveToFile (std::ofstream &file) const;
        template <class T, class... Args>
        static std::shared_ptr<T> &create (Args &&...args) {
            return std::make_shared<T>(std::forward<Args>(args)...);
        }
    protected:
        GUIType m_type;
        sf::Vector2f m_size;
        bool m_visible, m_resizable, m_needBorders;
        static Line m_line;
        static std::map<std::string, std::function<std::shared_ptr<GUIElement>()>> &getFabric ();
        template <class T>
        struct GuiFabric {
            GuiFabric (const std::string &str) {
                if (str.empty()) {
                    throw std::logic_error("GuiFabric::GuiFabric: empty key not allowed");
                }
                auto &fabric = GUIElement::getFabric();
                std::cout << "add " << str << " size: " << sizeof(T) << '\n';
                fabric.insert(std::make_pair(str, [] () -> std::shared_ptr<GUIElement> {
                    return std::make_shared<T>();
                }));
            }
        };
};

#endif