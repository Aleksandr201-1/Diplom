#ifndef SELECTOR_HPP
#define SELECTOR_HPP

#include <GUI/Interactable.hpp>
#include <GUI/TextBased.hpp>
#include <GUI/RenderTextureBased.hpp>

class Selector : public Interactable, public TextBased, public RenderTextureBased {
    private:
        void recalculateSize () override;
        void updateRender () override;
        void updateContent () override;
        void interact (const sf::Event &event, const sf::Vector2f &mousePos) override;
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
    public:
        Selector ();
        Selector (const std::vector<sf::String> &list);
        ~Selector ();
        void addVar (const sf::String &str);
        void setVarList (const std::vector<sf::String> &list);
        const sf::String &getCurrVal () const;
        uint64_t getCurrIdx () const;
        void addContentToSize (bool status);

        void setVisibleOffset (uint64_t offset);
        uint64_t getVisibleOffset () const;
        //void update () override;

        using GUIElement::setSize;
        void setSize (const sf::Vector2f &size) override;
        void loadFromFile (std::ifstream &file) override;
        void saveToFile (std::ofstream &file) const override;

        //void update (const sf::Event &event, const sf::Vector2f &mousePos);
    private:
        static const std::string M_NAME;
        static const GuiFabric<Selector> m_selectorFabric;
        static const uint64_t CONTENT_CAP;
        static const uint64_t MAX_LENGTH;
        static const uint64_t MIN_LENGTH;
        bool m_countContentSize;
        float m_borderOffset;
        uint64_t visibleOffset, m_currIdx, m_checkIdx;
        //sf::RectangleShape currentRect, visibleRect, checkRect;
        sf::RectangleShape checkRect;
        std::vector<sf::String> vars;
        sf::Sprite m_currentContent;
};

#endif