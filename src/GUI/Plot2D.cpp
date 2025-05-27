#include <GUI/Plot2D.hpp>

const std::string Plot2D::M_NAME = "Plot2D";
const Plot2D::GuiFabric<Plot2D> Plot2D::m_plot2DFabric(Plot2D::M_NAME);
const uint64_t Plot2D::M_ITERATION_LIMIT = 200;
const sf::Vector2f Plot2D::M_OFFSET(100.f, 20.f);

Plot2D::GraphInfo::GraphInfo () {}

Plot2D::GraphInfo::GraphInfo (const std::function<float (float)> &func, const Line &line, const sf::String &name) : func(func), line(line), name(name) {}

//Plot2D::GraphInfo::GraphInfo (const std::vector<sf::Vector2f> &points, const sf::Color &color) : color(color) {}

Plot2D::GraphInfo::~GraphInfo () {}

void Plot2D::recalculateSize () {

}

void Plot2D::updateRender () {
    std::ostringstream sstream;
    sf::Vertex lines[2];
    sstream.precision(2);
    float textWidth, textHeight;
    //auto renderSize = render.getSize();
    sf::RectangleShape rect;
    auto renderSize = sf::Vector2f(m_size.x - M_OFFSET.x * 2, m_size.y - M_OFFSET.x);
    rect.setSize(renderSize);
    rect.setFillColor(sf::Color::White);
    rect.setPosition(M_OFFSET);
    //clearing
    m_render.clear(m_clearColor);
    //отрисовка границ
    m_render.draw(rect);
    lines[0].color = lines[1].color = sf::Color::Black;
    lines[0].position = {0, 0};
    lines[1].position = {renderSize.x, 0};
    addOffset(lines, M_OFFSET);
    m_render.draw(lines, 2, sf::Lines);
    lines[0].position = {1, 0};
    lines[1].position = {0, renderSize.y};
    addOffset(lines, M_OFFSET);
    m_render.draw(lines, 2, sf::Lines);
    lines[0].position = {0, renderSize.y - 1};
    lines[1].position = {renderSize.x, renderSize.y - 1};
    addOffset(lines, M_OFFSET);
    m_render.draw(lines, 2, sf::Lines);
    lines[0].position = {renderSize.x, 0};
    lines[1].position = {renderSize.x, renderSize.y};
    addOffset(lines, M_OFFSET);
    m_render.draw(lines, 2, sf::Lines);
    //отрисовка функций
    if (m_mode & LOG_X) {
        coeffs[START_X] = coeffs[END_X] / std::pow(10.f, m_splitX);
        setAxisScales();
    }
    float stepX = (coeffs[END_X] - coeffs[START_X]) / M_ITERATION_LIMIT;
    float x1, x2;
    for (uint64_t i = 0; i < m_graphs.size(); ++i) {
        Line l = m_graphs[i].line;
        for (uint64_t j = 1; j < M_ITERATION_LIMIT; ++j) {
            if (m_mode & LOG_X) {
                float maxLog = std::log10(std::abs(coeffs[END_X]));
                float minLog = std::log10(std::abs(coeffs[START_X]));
                uint64_t pps = (M_ITERATION_LIMIT / m_splitX);
                uint64_t nSec = j / pps;
                float left = coeffs[END_X] / std::pow(10.f, m_splitX - nSec);
                float right = coeffs[END_X] / std::pow(10.f, m_splitX - nSec - 1);
                float step = (right - left) / pps;
                x1 = j == 1 ? left : x2;
                //x1 = left + (j - 1) % pps * step;
                x2 = left + j % pps * step;
                l.position[0] = {x1, m_graphs[i].func(x1)};
                l.position[1] = {x2, m_graphs[i].func(x2)};
                l.position[0].x = coeffs[START_X] + (std::log10(std::abs(l.position[0].x)) - minLog) / (maxLog - minLog) * (coeffs[END_X] - coeffs[START_X]);
                l.position[1].x = coeffs[START_X] + (std::log10(std::abs(l.position[1].x)) - minLog) / (maxLog - minLog) * (coeffs[END_X] - coeffs[START_X]);
                //std::cout << "1: " << x1 << " 2: " << x2 << "\n";
                //std::cout << j << ": 1: " << l.position[0].x << " 2: " << l.position[1].x << "\n";
                //std::cout << j << ": 1: " << (l.position[0].x - coeffs[START_X]) * scale << " 2: " << (l.position[1].x - coeffs[START_X]) * scale << "\n";
                //l.position[0].x = (l.position[0].x - coeffs[START_X]) * scale;
                //l.position[1].x = (l.position[1].x - coeffs[START_X]) * scale;
            } else {
                x1 = (j - 1) * stepX + coeffs[START_X];
                x2 = j * stepX + coeffs[START_X];
                l.position[0] = {x1, m_graphs[i].func(x1)};
                l.position[1] = {x2, m_graphs[i].func(x2)};
            }
            l.position[0].x = (l.position[0].x - coeffs[START_X]) * coeffs[SCALE_X];
            l.position[1].x = (l.position[1].x - coeffs[START_X]) * coeffs[SCALE_X];
            
            if (m_mode & LOG_Y) {
                float maxLog = std::log10(coeffs[END_Y] - coeffs[START_Y]);
                l.position[0].y *= -std::log10(l.position[0].y);
                l.position[1].y *= -std::log10(l.position[1].y);
            }
            l.position[0].y = (-l.position[0].y + coeffs[END_Y] - coeffs[START_Y]*0) * coeffs[SCALE_Y];
            l.position[1].y = (-l.position[1].y + coeffs[END_Y] - coeffs[START_Y]*0) * coeffs[SCALE_Y];
            l.position[0].y = std::max(std::min(l.position[0].y, renderSize.y), 0.f);
            l.position[1].y = std::max(std::min(l.position[1].y, renderSize.y), 0.f);
            l.position[0] += M_OFFSET;
            l.position[1] += M_OFFSET;
            m_render.draw(l);
        }
    }

    lines[0].color = lines[1].color = sf::Color::Black;

    sstream << coeffs[START_X];
    m_text.setString(sstream.str());
    sstream.str(std::string());
    textWidth = m_text.getLocalBounds().width;
    m_text.setPosition({-textWidth / 2, renderSize.y + 5});
    m_text.move(M_OFFSET);
    m_render.draw(m_text);
    sstream << std::scientific;
    long long int degree = 1000;
    for (uint64_t i = 0; i < m_splitX; ++i) {
        float splitStepX = (coeffs[END_X] - coeffs[START_X]) / m_splitX;
        float num = coeffs[START_X] + splitStepX * (i + 1);
        sstream << num;
        //std::cout << sstream.str().substr(sstream.str().find('e') + 1) << '\n';
        degree = std::min(std::stoll(sstream.str().substr(sstream.str().find('e') + 1)), degree);
        sstream.str(std::string());
    }
    sstream.str(std::string());
    //sstream << std::fixed;
    float divid = std::pow(0.1f, -degree);
    for (uint64_t i = 0; i < m_splitX; ++i) {
        float num;
        float sectorSize = renderSize.x / m_splitX;
        float splitStepX = (coeffs[END_X] - coeffs[START_X]) / m_splitX;
        if (m_mode & LOG_X) {
            num = coeffs[START_X] * std::pow(10.f, i + 1);
        } else {
            num = coeffs[START_X] + splitStepX * (i + 1);
        }
        sstream << num;
        m_text.setString(sstream.str());
        sstream.str(std::string());
        textWidth = m_text.getLocalBounds().width;
        m_text.setPosition({sectorSize * (i + 1) - textWidth / 2, renderSize.y + 5});
        m_text.move(M_OFFSET);
        for (uint64_t j = 1; j < 5; ++j) {
            lines[0].position = {sectorSize * i + (sectorSize / 5) * j, renderSize.y};
            lines[1].position = {sectorSize * i + (sectorSize / 5) * j, renderSize.y - 5};
            addOffset(lines, M_OFFSET);
            m_render.draw(lines, 2, sf::Lines);
        }
        lines[0].position = {sectorSize * (i + 1), renderSize.y};
        lines[1].position = {sectorSize * (i + 1), renderSize.y - 10};
        addOffset(lines, M_OFFSET);
        m_render.draw(lines, 2, sf::Lines);
        m_render.draw(m_text);
    }
    m_text.setString(m_XName);
    textWidth = m_text.getLocalBounds().width;
    m_text.setPosition(sf::Vector2f(renderSize.x / 2 - textWidth / 2, renderSize.y + 40) + M_OFFSET);
    m_render.draw(m_text);

    sstream << coeffs[END_Y];
    m_text.setString(sstream.str());
    sstream.str(std::string());
    textWidth = m_text.getLocalBounds().width;
    textHeight = m_text.getLocalBounds().height;
    m_text.setPosition({-textWidth - 10, -15});
    m_text.move(M_OFFSET);
    m_render.draw(m_text);
    for (uint64_t i = 0; i < m_splitY; ++i) {
        float sectorSize = renderSize.y / m_splitY;
        float splitStepY = (coeffs[END_Y] - coeffs[START_Y]) / m_splitY;
        float num = coeffs[END_Y] - splitStepY * (i + 1);
        sstream << num;
        m_text.setString(sstream.str());
        sstream.str(std::string());
        textWidth = m_text.getLocalBounds().width;
        m_text.setPosition({-textWidth - 10, sectorSize * (i + 1) - 15});
        m_text.move(M_OFFSET);
        for (uint64_t j = 1; j < 5; ++j) {
            lines[0].position = {0, sectorSize * i + (sectorSize / 5) * j};
            lines[1].position = {5, sectorSize * i + (sectorSize / 5) * j};
            addOffset(lines, M_OFFSET);
            m_render.draw(lines, 2, sf::Lines);
        }
        lines[0].position = {0, sectorSize * (i + 1)};
        lines[1].position = {10, sectorSize * (i + 1)};
        addOffset(lines, M_OFFSET);
        m_render.draw(lines, 2, sf::Lines);
        m_render.draw(m_text);
    }
    m_text.setString(m_YName);
    textWidth = m_text.getLocalBounds().width;
    m_text.rotate(-90.f);
    m_text.setPosition(sf::Vector2f(0.f, renderSize.y / 2 + textWidth / 2));
    m_render.draw(m_text);
    m_text.setRotation(0.f);

    //отрисовка сетки
    if (m_gridX.has_value() && !(m_mode & ZONES)) {
        Line gridLine(*m_gridX, sf::Color(127, 127, 127, 127), 1.f, 10.f);
        for (uint64_t i = 0; i < m_splitX - 1; ++i) {
            float sectorSize = renderSize.x / m_splitX;
            gridLine.position[0] = {sectorSize * (i + 1), renderSize.y};
            gridLine.position[1] = {sectorSize * (i + 1), 0};
            gridLine.position[0] += M_OFFSET;
            gridLine.position[1] += M_OFFSET;
            m_render.draw(gridLine);
        }
    }
    if (m_gridY.has_value() && !(m_mode & ZONES)) {
        Line gridLine(*m_gridX, sf::Color(127, 127, 127, 127), 1.f, 10.f);
        for (uint64_t i = 0; i < m_splitY - 1; ++i) {
            float sectorSize = renderSize.y / m_splitY;
            gridLine.position[0] = {0, sectorSize * (i + 1)};
            gridLine.position[1] = {renderSize.x, sectorSize * (i + 1)};
            gridLine.position[0] += M_OFFSET;
            gridLine.position[1] += M_OFFSET;
            m_render.draw(gridLine);
        }
    }

    m_render.display();
    m_needRedraw = false;
}

void Plot2D::updateContent () {
    if (m_needRedraw) {
        updateRender();
    }
    m_needUpdate = false;
    //m_restart.update();
}

void Plot2D::interact (const sf::Event &event, const sf::Vector2f &mousePos) {
    auto newMousePos = mousePos - M_OFFSET;
    auto renderSize = sf::Vector2f(m_size.x - M_OFFSET.x * 2, m_size.y - M_OFFSET.x);
    if (newMousePos.x > 0 && newMousePos.x < renderSize.x && newMousePos.y > 0 && newMousePos.y < renderSize.y) {
        sf::Vector2f currPos;
        currPos.x = newMousePos.x / coeffs[SCALE_X] + coeffs[START_X];
        currPos.y = -newMousePos.y / coeffs[SCALE_Y] + coeffs[END_Y];
        if (!m_graphs.empty()) {
            uint64_t idx = 0;
            float currY = m_graphs[idx].func(currPos.x);
            float diff = std::abs(currY - currPos.y);
            for (uint64_t i = 1; i < m_graphs.size(); ++i) {
                float tmp = m_graphs[i].func(currPos.x);
                float currDiff = std::abs(tmp - currPos.y);
                if (currDiff < diff) {
                    diff = currDiff;
                    currY = tmp;
                }
            }
            std::ostringstream ss;
            ss.precision(4);
            ss << std::scientific;
            ss << "(" << currPos.x << "; " << currY << ")";
            m_pointCoords = sf::String(ss.str());
            // m_circle.setRadius(2.f);
            // m_circle.setFillColor(sf::Color::White);
            // m_circle.setOutlineColor(sf::Color::Black);
            // m_circle.setOutlineThickness(1.f);
            // m_circle.setOrigin({1.f, 1.f});
            // m_circle.setPosition(sf::Vector2f{newMousePos.x, (-currY + coeffs[END_Y]) * coeffs[SCALE_Y]} + M_OFFSET);
            m_text.setString(m_pointCoords);
            //m_text.setPosition(m_circle.getPosition() + sf::Vector2f(5, 5));
        }
        //m_text.setString("(" + toString(currPos.x, 2) + "; " + toString(currPos.y, 2) + ")");
        if (event.type == sf::Event::MouseWheelMoved) {
            coeffs[START_X] += event.mouseWheel.delta * 5 / coeffs[SCALE_X];
            coeffs[START_Y] += event.mouseWheel.delta * 5 / coeffs[SCALE_Y];
            coeffs[END_X] -= event.mouseWheel.delta * 5 / coeffs[SCALE_X];
            coeffs[END_Y] -= event.mouseWheel.delta * 5 / coeffs[SCALE_Y];
            setAxisScales();
            m_needUpdate = m_needRedraw = true;
        }
        if (m_interactFocus) {
            if (event.type == sf::Event::MouseMoved && sf::Mouse::isButtonPressed(sf::Mouse::Left)) {
                float distX = coeffs[END_X] - coeffs[START_X], distY = coeffs[END_Y] - coeffs[START_Y];
                coeffs[START_X] += (m_oldMousePos.x - newMousePos.x) / coeffs[SCALE_X];
                coeffs[END_X] += (m_oldMousePos.x - newMousePos.x) / coeffs[SCALE_X];
                coeffs[START_Y] += (newMousePos.y - m_oldMousePos.y) / coeffs[SCALE_Y];
                coeffs[END_Y] += (newMousePos.y - m_oldMousePos.y) / coeffs[SCALE_Y];
                m_needUpdate = m_needRedraw = true;
            }
        }
    }
    m_oldMousePos = newMousePos;
}

void Plot2D::drawVisible (sf::RenderTarget& target, sf::RenderStates states) const {
    target.draw(m_visibleContent, states);
    //target.draw(m_circle, states);
    target.draw(m_text, states);
}

sf::Color Plot2D::getRandomColor () {
    uint32_t num = static_cast<uint32_t>(std::rand());
    num = (num << 8) + 255;
    return sf::Color(num);
}

// void Plot2D::setAxisBorders (uint64_t start) {
//     // for (uint64_t i = start; i < m_graphs.size(); ++i) {
//     //     for (uint64_t j = 0; j < M_ITERATION_LIMIT; ++j) {
//     //         coeffs[OLD_START_X] = coeffs[START_X] = std::min(coeffs[START_X], m_graphs[i].points[j].x);
//     //         coeffs[OLD_END_X] = coeffs[END_X] = std::max(m_graphs[i].points[j].x, coeffs[END_X]);
//     //         coeffs[OLD_START_Y] = coeffs[START_Y] = std::min(coeffs[START_Y], m_graphs[i].points[j].y);
//     //         coeffs[OLD_END_Y] = coeffs[END_Y] = std::max(m_graphs[i].points[j].y, coeffs[END_Y]);
//     //     }
//     // }
// }

void Plot2D::setAxisBorders (const std::vector<sf::Vector2f> &points) {
    for (uint64_t i = 0; i < points.size(); ++i) {
        coeffs[OLD_START_X] = coeffs[START_X] = std::min(coeffs[START_X], points[i].x);
        coeffs[OLD_END_X] = coeffs[END_X] = std::max(points[i].x, coeffs[END_X]);
        coeffs[OLD_START_Y] = coeffs[START_Y] = std::min(coeffs[START_Y], points[i].y);
        coeffs[OLD_END_Y] = coeffs[END_Y] = std::max(points[i].y, coeffs[END_Y]);
    }
    if (coeffs[START_Y] == coeffs[END_Y]) {
        double y = coeffs[START_Y];
        coeffs[OLD_START_Y] = coeffs[START_Y] -= std::abs(coeffs[START_Y]) * 2.5;
        coeffs[OLD_END_Y] = coeffs[END_Y] += std::abs(coeffs[START_Y]) * 2.5;
    }
}

void Plot2D::setAxisBorders (const std::vector<float> &Xpoints, const std::vector<float> &Ypoints) {
    for (uint64_t i = 0; i < Xpoints.size(); ++i) {
        coeffs[OLD_START_X] = coeffs[START_X] = std::min(coeffs[START_X], Xpoints[i]);
        coeffs[OLD_END_X] = coeffs[END_X] = std::max(Xpoints[i], coeffs[END_X]);
        coeffs[OLD_START_Y] = coeffs[START_Y] = std::min(coeffs[START_Y], Ypoints[i]);
        coeffs[OLD_END_Y] = coeffs[END_Y] = std::max(Ypoints[i], coeffs[END_Y]);
    }
    if (coeffs[START_Y] == coeffs[END_Y]) {
        double y = coeffs[START_Y];
        coeffs[OLD_START_Y] = coeffs[START_Y] -= std::abs(coeffs[START_Y]) * 2.5;
        coeffs[OLD_END_Y] = coeffs[END_Y] += std::abs(coeffs[START_Y]) * 2.5;
    }
}

void Plot2D::setAxisScales () {
    //scaleX = ((float)render.getSize().x) / (endX - startX);
    //scaleY = ((float)render.getSize().y) / (endY - startY);
    auto renderSize = sf::Vector2f(m_size.x - M_OFFSET.x * 2, m_size.y - M_OFFSET.x);
    coeffs[OLD_SCALE_X] = coeffs[SCALE_X] = renderSize.x / (coeffs[END_X] - coeffs[START_X]);
    coeffs[OLD_SCALE_Y] = coeffs[SCALE_Y] = renderSize.y / (coeffs[END_Y] - coeffs[START_Y]);
}

void Plot2D::addOffset (sf::Vertex lines[2], const sf::Vector2f &off) {
    lines[0].position += off;
    lines[1].position += off;
}

Plot2D::Plot2D () {
    m_needUpdate = m_needRedraw = true;
    m_needTriggerPressed = true;
    m_pointInfo = true;
    m_oldMousePos = {0, 0};
    m_size = {600, 500};
    setRenderSize(600, 500);
    m_text.setCharacterSize(15);
    m_gridX = m_gridY = std::make_optional<Line::Type>(Line::Type::DOTTED);
    m_XName = "X";
    m_YName = "Y";
    m_splitX = m_splitY = 10;
    m_mode = FUNCTIONAL;
    m_clearColor = sf::Color::Transparent;
    // coeffs[OLD_START_X] = coeffs[START_X] = 0;
    // coeffs[OLD_END_X] = coeffs[END_X] = 10;
    // coeffs[OLD_START_Y] = coeffs[START_Y] = 0;
    // coeffs[OLD_END_Y] = coeffs[END_Y] = 10;
    coeffs[START_X] = coeffs[START_Y] = INFINITY;
    coeffs[END_X] = coeffs[END_Y] = -INFINITY;
    //setAxisScales();
}

Plot2D::Plot2D (const std::vector<std::function<float (float)>> &funcs) : Plot2D() {
    std::vector<Line> lines(funcs.size());
    std::vector<sf::String> names(funcs.size());
    for (uint64_t i = 0; i < funcs.size(); ++i) {
        lines[i].type = Line::Type(i % 4);
        lines[i].color = getRandomColor();
        lines[i].thickness = 2.5f;
        names[i] = sf::String(std::to_string(i + 1));
    }
    setFuncs(funcs, lines, names);
}

Plot2D::Plot2D (const std::vector<std::vector<sf::Vector2f>> &points) : Plot2D() {
    std::vector<Line> lines(points.size());
    std::vector<sf::String> names(points.size());
    for (uint64_t i = 0; i < points.size(); ++i) {
        lines[i].type = Line::Type(i % 4);
        lines[i].color = getRandomColor();
        lines[i].thickness = 2.5f;
        names[i] = sf::String(std::to_string(i + 1));
    }
    setFuncs(points, lines, names);
}

Plot2D::Plot2D (const std::vector<std::function<float (float)>> &funcs, const std::vector<Line> &lines, const std::vector<sf::String> &names) : Plot2D() {
    setFuncs(funcs, lines, names);
}

Plot2D::Plot2D (const std::vector<std::vector<sf::Vector2f>> &points, const std::vector<Line> &lines, const std::vector<sf::String> &names) : Plot2D() {
    setFuncs(points, lines, names);
}

Plot2D::~Plot2D () {}

void Plot2D::addFunc (const std::function<float (float)> &func, const Line &line, const sf::String &name) {
    m_needUpdate = m_needRedraw = true;
    // float step = (x2 - x1) / iterLimit;
    // std::vector<sf::Vector2f> points(iterLimit);
    // for (uint64_t i = 0; i < iterLimit; ++i) {
    //     points[i].x = x1 + step * i;
    //     points[i].y = func(points[i].x);
    // }
    m_graphs.push_back(GraphInfo(func, line, name));
    //setAxisBorders(points);
}

void Plot2D::addFunc (const std::vector<sf::Vector2f> &points, const Line &line, const sf::String &name) {
    m_needUpdate = m_needRedraw = true;
    auto func = [=] (float x) -> float {
        return LinearInterpolation(points, x);
    };
    addFunc(func, line, name);
    setAxisBorders(points);
    setAxisScales();
}

void Plot2D::addFunc (const std::vector<float> &Xpoints, const std::vector<float> &Ypoints, const Line &line, const sf::String &name) {
    m_needUpdate = m_needRedraw = true;
    auto func = [=] (float x) -> float {
        return LinearInterpolation(Xpoints, Ypoints, x);
    };
    addFunc(func, line, name);
    setAxisBorders(Xpoints, Ypoints);
    setAxisScales();
}

void Plot2D::setFuncs (const std::vector<std::function<float (float)>> &funcs, const std::vector<Line> &lines, const std::vector<sf::String> &names) {
    m_needUpdate = m_needRedraw = true;
    float stepX = (coeffs[END_X] - coeffs[START_X]) / M_ITERATION_LIMIT;
    m_graphs.resize(funcs.size());
    for (uint64_t i = 0; i < funcs.size(); ++i) {
        m_graphs[i].func = funcs[i];
        for (uint64_t j = 0; j <= M_ITERATION_LIMIT; ++j) {
            float currY =  funcs[i](coeffs[START_X] + j * stepX);
            coeffs[OLD_START_Y] = coeffs[START_Y] = std::min(coeffs[START_Y], currY);
            coeffs[OLD_END_Y] = coeffs[END_Y] = std::max(currY, coeffs[END_Y]);
        }
        m_graphs[i].line = lines[i];
        m_graphs[i].name = names[i];
    }
    setAxisScales();
}

void Plot2D::setFuncs (const std::vector<std::vector<sf::Vector2f>> &points, const std::vector<Line> &lines, const std::vector<sf::String> &names) {
    m_needUpdate = m_needRedraw = true;
    clear();
    m_graphs.resize(points.size());
    //m_polynoms.resize(points.size());
    //m_points.resize(points.size());
    for (uint64_t i = 0; i < points.size(); ++i) {
        //m_polynoms[i] = CubeSpline(points[i]);
        //m_points[i] = points[i];
        m_graphs[i].func = [=] (float x) -> float {
            return LinearInterpolation(points[i], x);
        };
        m_graphs[i].line = lines[i];
        m_graphs[i].name = names[i];
        setAxisBorders(points[i]);
    }
    setAxisScales();
}

void Plot2D::setXBorders (float x1, float x2) {
    coeffs[OLD_START_X] = coeffs[START_X] = x1;
    coeffs[OLD_END_X] = coeffs[END_X] = x2;
}

void Plot2D::setYBorders (float y1, float y2) {
    coeffs[OLD_START_Y] = coeffs[START_Y] = y1;
    coeffs[OLD_END_Y] = coeffs[END_Y] = y2;
}

void Plot2D::setMap (const std::vector<std::vector<float>> &map) {
    //m_polynoms = map;
}

void Plot2D::setMode (Mode mode) {
    m_mode = mode;
}

void Plot2D::addMode (Mode mode) {
    m_mode |= uint64_t(mode);
}

void Plot2D::setColor (uint64_t i, const sf::Color &color) {
    m_needUpdate = m_needRedraw = true;
    m_graphs[i].line.color = color;
}

void Plot2D::setColors (const std::vector<sf::Color> &colors) {
    m_needUpdate = m_needRedraw = true;
    uint64_t min = std::min(colors.size(), m_graphs.size());
    for (uint64_t i = 0; i < min; ++i) {
        m_graphs[i].line.color = colors[i];
    }
}

void Plot2D::setName (uint64_t i, const sf::String &name) {
    m_needUpdate = m_needRedraw = true;
    m_graphs[i].name = name;
}

void Plot2D::setNames (const std::vector<sf::String> &names) {
    m_needUpdate = m_needRedraw = true;
    uint64_t min = std::min(names.size(), m_graphs.size());
    for (uint64_t i = 0; i < min; ++i) {
        m_graphs[i].name = names[i];
    }
}

void Plot2D::setLine (uint64_t i, const Line &line) {
    m_needUpdate = m_needRedraw = true;
    m_graphs[i].line = line;
}

void Plot2D::setLines (const std::vector<Line> &lines) {
    m_needUpdate = m_needRedraw = true;
    uint64_t min = std::min(lines.size(), m_graphs.size());
    for (uint64_t i = 0; i < min; ++i) {
        m_graphs[i].line = lines[i];
    }
}

void Plot2D::setXName (const sf::String &XName) {
    m_needUpdate = m_needRedraw = true;
    m_XName = XName;
}

void Plot2D::setYName (const sf::String &YName) {
    m_needUpdate = m_needRedraw = true;
    m_YName = YName;
}

void Plot2D::needPointInfo (bool status) {
    m_pointInfo = status;
}

sf::Vector2f Plot2D::getCursorPos () const {
    return m_oldMousePos;
}

std::shared_ptr<Legend> Plot2D::generateLegend () {
    auto legend = std::make_shared<Legend>();
    for (uint64_t i = 0; i < m_graphs.size(); ++i) {
        Legend::Record rec(m_graphs[i].line, m_graphs[i].name);
        legend->add(rec);
    }
    return legend;
}

void Plot2D::clear () {
    coeffs[START_X] = coeffs[START_Y] = INFINITY;
    coeffs[END_X] = coeffs[END_Y] = -INFINITY;
    m_graphs.clear();
    //m_polynoms.clear();
}

void Plot2D::setSize (const sf::Vector2f &size) {
    if (size.x < 200 || size.y < 100) {
        return;
    }
    m_needUpdate = m_needRedraw = true;
    m_oldMousePos = {0, 0};
    m_size = size;
    setRenderSize(m_size.x, m_size.y);
    setAxisScales();
}

void Plot2D::loadFromFile (std::ifstream &file) {
    Interactable::loadFromFile(file);
    TextBased::loadFromFile(file);
    RenderTextureBased::loadFromFile(file);
}

void Plot2D::saveToFile (std::ofstream &file) const {
    uint64_t size = M_NAME.size();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    file.write(&M_NAME[0], sizeof(char) * size);
    file.write(reinterpret_cast<const char*>(coeffs), sizeof(coeffs[0]) * AXIS_COEFFS::COUNT);
    file.write(reinterpret_cast<const char*>(m_mode), sizeof(m_mode));
    file.write(reinterpret_cast<const char*>(&m_splitX), sizeof(m_splitX));
    file.write(reinterpret_cast<const char*>(&m_splitY), sizeof(m_splitY));
    //file.write(reinterpret_cast<const char*>(&m_needGridX), sizeof(m_needGridX));
    //file.write(reinterpret_cast<const char*>(&m_needGridY), sizeof(m_needGridY));
    file.write(reinterpret_cast<const char*>(&m_lockBorders), sizeof(m_lockBorders));
    size = m_XName.getSize();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    if (size > 0) {
        file.write(reinterpret_cast<const char*>(m_XName.getData()), sizeof(m_XName[0]) * size);
    }
    size = m_YName.getSize();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    if (size > 0) {
        file.write(reinterpret_cast<const char*>(m_YName.getData()), sizeof(m_YName[0]) * size);
    }
    file.write(reinterpret_cast<const char*>(&m_plotSize), sizeof(m_plotSize));
    // file.write(reinterpret_cast<const char*>(&status), sizeof(status));
    // file.write(reinterpret_cast<const char*>(&status), sizeof(status));
    // file.write(reinterpret_cast<const char*>(&status), sizeof(status));
    // float coeffs[START_X], coeffs[END_X], coeffs[START_Y], coeffs[END_Y], stepX, m_stepY, coeffs[SCALE_X], coeffs[SCALE_Y];
    // float coeffs[OLD_START_X], coeffs[OLD_START_Y], coeffs[OLD_END_X], coeffs[OLD_END_Y];
    // float coeffs[OLD_SCALE_X], coeffs[OLD_SCALE_Y];
    // uint64_t m_mode, m_splitX, m_splitY;
    // bool m_needGrid, m_needLegend, m_showPoint;


    Interactable::saveToFile(file);
    TextBased::saveToFile(file);
    RenderTextureBased::saveToFile(file);
}