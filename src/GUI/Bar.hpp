#ifndef BAR_HPP
#define BAR_HPP

#include <GUI/Interactable.hpp>

class Bar : public Interactable {
    private:
        void updateContent () override;
        void interact (const sf::Event &event, const sf::Vector2f &mousePos) override;
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
    public:
        Bar ();
        //Bar (uint64_t min = 0, uint64_t max = UINT64_MAX);
        Bar (uint64_t &param, uint64_t min = 0, uint64_t max = UINT64_MAX);
        ~Bar ();
        //void setSprites (const SpriteManager &manager, uint64_t first, uint64_t second);
        void setSprites (const sf::Sprite &box, const sf::Sprite &bar);
        void setRef(uint64_t &param);
        void setMinMax (uint64_t min = 0, uint64_t max = UINT64_MAX);
        using GUIElement::setSize;
        void setSize (const sf::Vector2f &size) override;
        void loadFromFile (std::ifstream &file) override;
        void saveToFile (std::ofstream &file) const override;
    private:
        static const std::string M_NAME;
        static const GuiFabric<Bar> M_BAR_FABRIC;
        std::reference_wrapper<uint64_t> m_param;
        //const uint64_t& m_param;
        uint64_t m_min, m_max;
        sf::Sprite m_box, m_bar;
        //sf::RectangleShape m_box, m_bar;
        sf::IntRect m_rect;
};

#endif