#ifndef SLIDER_HPP
#define SLIDER_HPP

#include <GUI/Interactable.hpp>
//#include <General/SpriteManager.hpp>
#include <cmath>

class Slider : public Interactable {
    private:
        void interact (const sf::Event &event, const sf::Vector2f &mousePos) override;
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
    public:
        //Slider (const SpriteManager &manager, uint64_t first, uint64_t second, uint64_t min, uint64_t max);
        Slider ();
        Slider (float min, float max, uint64_t division = 1);
        ~Slider ();
        void setValue (float value);
        float getValue () const;
        void setMinMax (float min, float max);
        void setDivision (uint64_t division);
        void setSliderSize (float size);
        void setSensativity (float sensativity);
        void roundValue (bool status);
        using GUIElement::setSize;
        void setSize (const sf::Vector2f &size) override;
        void loadFromFile (std::ifstream &file) override;
        void saveToFile (std::ofstream &file) const override;
    private:
        static const std::string M_NAME;
        static const GuiFabric<Slider> m_sliderFabric;
        bool m_needRound;
        uint64_t m_division;
        float m_value, m_sensativity, m_minimalStep;
        float m_min, m_max;
        sf::RectangleShape m_box, m_slider;
};

#endif