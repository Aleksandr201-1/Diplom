#ifndef PLOT2D_HPP
#define PLOT2D_HPP

#include <General/FloatToString.hpp>
#include <NumericMethods/Interpolation.hpp>
#include <SFExtensions/Line.hpp>
#include <GUI/Interactable.hpp>
#include <GUI/RenderTextureBased.hpp>
#include <GUI/TextBased.hpp>
#include <GUI/Legend.hpp>
#include <functional>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <optional>

//класс для отрисовки графика на двухмерной плоскости
class Plot2D : public Interactable, public TextBased, public RenderTextureBased, public LegendAvailable {
    private:
        enum AXIS_COEFFS {
            START_X = 0,
            END_X,
            SCALE_X,
            OLD_START_X,
            OLD_END_X,
            OLD_SCALE_X,
            START_Y,
            END_Y,
            SCALE_Y,
            OLD_START_Y,
            OLD_END_Y,
            OLD_SCALE_Y,
            COUNT
        };
    public:
        enum Mode : uint64_t {
            FUNCTIONAL  = 1 << 0,
            POINTS      = 1 << 1,
            PARAMETRIC  = 1 << 2,
            LOGARITHMIC = 1 << 3,
            ZONES       = 1 << 4,
            LOG_X       = 1 << 5,
            LOG_Y       = 1 << 6
        };
    private:
        struct GraphInfo {
            GraphInfo ();
            GraphInfo (const GraphInfo &) = default;
            GraphInfo (GraphInfo &&) = default;
            GraphInfo (const std::function<float (float)> &func, const Line &line, const sf::String &name = L"");
            ~GraphInfo ();

            std::function<float (float)> func;
            Line line;
            sf::String name;
        };
    private:
        void recalculateSize () override;
        void updateRender () override;
        void updateContent () override;
        void interact (const sf::Event &event, const sf::Vector2f &mousePos) override;
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
        sf::Color getRandomColor ();
        void setAxisBorders (const std::vector<sf::Vector2f> &points);
        void setAxisBorders (const std::vector<float> &Xpoints, const std::vector<float> &Ypoints);
        void setAxisScales ();
        void addOffset (sf::Vertex lines[2], const sf::Vector2f &off);
    public:
        Plot2D ();
        Plot2D (const Plot2D &plot) = default;
        Plot2D (Plot2D &&plot) = default;
        Plot2D (const std::vector<std::function<float (float)>> &funcs);
        Plot2D (const std::vector<std::vector<sf::Vector2f>> &points);
        Plot2D (const std::vector<std::function<float (float)>> &funcs, const std::vector<Line> &lines, const std::vector<sf::String> &names);
        Plot2D (const std::vector<std::vector<sf::Vector2f>> &points, const std::vector<Line> &lines, const std::vector<sf::String> &names);
        ~Plot2D ();

        void addFunc (const std::function<float (float)> &func, const Line &line = Line(), const sf::String &name = L"");
        void addFunc (const std::vector<sf::Vector2f> &points, const Line &line = Line(), const sf::String &name = L"");
        void addFunc (const std::vector<float> &Xpoints, const std::vector<float> &Ypoints, const Line &line = Line(), const sf::String &name = L"");
        void setFuncs (const std::vector<std::function<float (float)>> &funcs, const std::vector<Line> &lines, const std::vector<sf::String> &names);
        void setFuncs (const std::vector<std::vector<sf::Vector2f>> &points, const std::vector<Line> &lines, const std::vector<sf::String> &names);
        void setXBorders (float x1, float x2);
        void setYBorders (float y1, float y2);
        void setMap (const std::vector<std::vector<float>> &map);
        void setMode (Mode mode);
        void addMode (Mode mode);
        void setColor (uint64_t i, const sf::Color &color);
        void setColors (const std::vector<sf::Color> &colors);
        void setName (uint64_t i, const sf::String &name);
        void setNames (const std::vector<sf::String> &names);
        void setLine (uint64_t i, const Line &line);
        void setLines (const std::vector<Line> &lines);
        void setXName (const sf::String &XName);
        void setYName (const sf::String &YName);
        void needPointInfo (bool status);
        sf::Vector2f getCursorPos () const;
        std::shared_ptr<Legend> generateLegend () override;
        void clear ();
        Plot2D &operator= (const Plot2D &plot) = default;
        //оператор перемещения
        Plot2D &operator= (Plot2D &&plot) = default;
        using GUIElement::setSize;
        void setSize (const sf::Vector2f &size) override;
        void loadFromFile (std::ifstream &file) override;
        void saveToFile (std::ofstream &file) const override;
    private:
        static const std::string M_NAME;
        static const GuiFabric<Plot2D> m_plot2DFabric;
        static const uint64_t M_ITERATION_LIMIT;
        static const sf::Vector2f M_OFFSET;
        static sf::CircleShape m_circle;
        static sf::RectangleShape m_rectangle;
        float coeffs[AXIS_COEFFS::COUNT];
        uint64_t m_mode, m_splitX, m_splitY;
        bool m_lockBorders, m_pointInfo;
        std::optional<Line::Type> m_gridX, m_gridY;
        sf::String m_XName, m_YName, m_pointCoords;
        std::vector<std::vector<sf::Vector2f>> m_points;
        sf::Vector2f m_plotSize;
        std::vector<GraphInfo> m_graphs;
        sf::Vector2f m_oldMousePos;
};

#endif