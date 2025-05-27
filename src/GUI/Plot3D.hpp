#ifndef PLOT3D_HPP
#define PLOT3D_HPP

#include <SFML/OpenGL.hpp>
#include <GUI/Interactable.hpp>
#include <GUI/TextBased.hpp>
#include <GUI/RenderTextureBased.hpp>
#include <General/FloatToString.hpp>
#include <functional>
#include <cmath>

enum class Plot3DMode : uint64_t {
    FUNCTIONAL,
    POINTS,
    PARAMETRIC,
    LOGARITHMIC
};

//класс для отрисовки графика
class Plot3D : public Interactable, public TextBased, public RenderTextureBased {
    private:
        void createNormals ();
        //void checkForForbiddenPoints ();
        void recalculateSize () override;
        void updateRender () override;
        void updateContent () override;
        void interact (const sf::Event &event, const sf::Vector2f &mousePos) override;
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
        sf::Color getRandomColor ();
        void setAxisBorders ();
        void setAxisScales ();
        //void addOffset (const sf::Vector2f &off);
    public:
        Plot3D ();
        Plot3D (const std::function<float (float, float)> &func, uint64_t iterLimit, float x1, float x2, float z1, float z2);
        Plot3D (const std::vector<sf::Vector3f> &points, uint64_t iterX, uint64_t iterZ);
        ~Plot3D ();

        void setPlane (const std::function<float (float, float)> &func, uint64_t iterLimit, float x1, float x2, float z1, float z2);
        void setPlane (const std::vector<sf::Vector3f> &points, uint64_t iterX, uint64_t iterZ);
        void setXBorders (float x1, float x2);
        void setYBorders (float y1, float y2);
        void setZBorders (float z1, float z2);
        void setLineColor (const sf::Color &color);
        void setPlaneColor (const sf::Color &color);
        void clear ();
        //void update () override;
        using GUIElement::setSize;
        void setSize (const sf::Vector2f &v) override;
    private:
        static const uint64_t ITERATION_LIMIT = 200;
        float startX, endX, scaleX, startY, endY, scaleY, startZ, endZ, scaleZ;
        float rotateX, rotateY, m_scale;
        uint64_t mode, splitX, splitY, splitZ;
        bool m_needGrid, m_needLegend, m_showPoints, m_needLight;
        sf::Color m_lineColor, m_planeColor;
        sf::Vector2f plotSize, currPos, oldMousePos;
        //std::vector<float> forbidden;
        //std::vector<std::vector<sf::Vector2f>> points;
        std::vector<sf::Vector3f> plane;
        std::vector<sf::Vector3f> normals;
};

#endif