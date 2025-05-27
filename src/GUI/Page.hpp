#ifndef PAGE_HPP
#define PAGE_HPP

#include <GUI/Interactable.hpp>
#include <GUI/GUIContainer.hpp>
#include <memory>
#include <map>

class Page : public GUIContainer {
    private:
        void recalculateSize () override;
    //     void updateContent () override;
    //     void interact (const sf::Event &event, const sf::Vector2f &mousePos) override;
    //     void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
    public:
        Page ();
        ~Page ();

        // void addGUIElement (uint64_t id, std::shared_ptr<GUIElement> &&el);
        // std::shared_ptr<GUIElement> &getGUIElement (uint64_t id);
        using GUIElement::setSize;
        void setSize (const sf::Vector2f &size) override;
        //void saveToFile (std::ofstream &file) const override;

        //std::shared_ptr<GUIElement> &operator[] (uint64_t id);
        //const std::shared_ptr<GUIElement> &operator[] (uint64_t id) const;
    //private:
        //static const std::string M_NAME;
        //std::map<uint64_t, std::shared_ptr<GUIElement>> m_content;
};

#endif