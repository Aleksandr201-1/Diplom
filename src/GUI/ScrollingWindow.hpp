#ifndef SCROLLING_WINDOW_HPP
#define SCROLLING_WINDOW_HPP

#include <GUI/Interactable.hpp>
#include <GUI/GUIContainer.hpp>
#include <GUI/RenderTextureBased.hpp>
#include <GUI/Slider.hpp>
#include <memory>

class ScrollingWindow : public GUIContainer, public RenderTextureBased {
    private:
        void recalculateSize () override;
        void updateRender () override;
        void updateContent () override;
        void interact (const sf::Event &event, const sf::Vector2f &mousePos) override;
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
    public:
        ScrollingWindow ();
        ~ScrollingWindow ();

        void setSensativity (float sensativity);
        void setWindowSize (uint64_t width, uint64_t height);
        void setWindowSize (const sf::Vector2u &size);
        void needXAxis (bool status);
        void needYAxis (bool status);
        //void update () override;
        using GUIElement::setSize;
        void setSize (const sf::Vector2f &size) override;
        void loadFromFile (std::ifstream &file) override;
        void saveToFile (std::ofstream &file) const override;
    private:
        static const std::string M_NAME;
        static const GuiFabric<ScrollingWindow> m_scrollingWindowFabric;
        static const uint64_t M_SLIDER_FREE_SPACE;
        Slider m_sliderX, m_sliderY;
        bool m_xAxis, m_yAxis;
        sf::Vector2f m_availableOffset;
        //std::map<uint64_t, std::shared_ptr<GUIElement>> m_gui;
};

#endif 