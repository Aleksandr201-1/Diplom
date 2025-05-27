#ifndef HEAT_MAP_HPP
#define HEAT_MAP_HPP

#include <SFML/OpenGL.hpp>
#include <GUI/Interactable.hpp>
#include <GUI/RenderTextureBased.hpp>
#include <Math/Matrix.hpp>
#include <cmath>

class HeatMap : public Interactable, public RenderTextureBased {
    private:
        void recalculateSize () override;
        void updateRender () override;
        void updateContent () override;
        void interact (const sf::Event &event, const sf::Vector2f &mousePos) override;
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
    public:
        HeatMap ();
        ~HeatMap ();

        void setMap (const Matrix<float> &map);

        using GUIElement::setSize;
        void setSize (const sf::Vector2f &size) override;
        void loadFromFile (std::ifstream &file) override;
        void saveToFile (std::ofstream &file) const override;
    private:
        Matrix<float> m_map;
        float m_min, m_max;
        float m_squareSizeX, m_squareSizeY;
};

#endif