#include <SFExtensions/Line.hpp>

void Line::draw (sf::RenderTarget& target, sf::RenderStates states) const {
    target.setActive(true);
    target.pushGLStates();
    states.transform *= getTransform();
    glViewport(0, 0, target.getSize().x, target.getSize().y);
    glLineWidth(thickness);
    glPointSize(thickness * 2);
    glLoadMatrixf(states.transform.getMatrix());
    sf::Vector2f thr1, thr2;
    uint64_t intervalCount = 0;
    float length = 0.f, xLen = position[0].x - position[1].x, yLen = position[0].y - position[1].y;
    length = std::sqrt(xLen * xLen + yLen * yLen);
    intervalCount = (uint64_t)(length / interval);
    switch (type) {
        case Type::POINT:
            if (interval == 0.f) {
                glBegin(GL_POINTS);
                    glColor3f((float) color.r / 255, (float) color.g / 255, (float) color.b / 255);
                    glVertex2f((position[0].x + position[1].x) / 2, (position[0].y + position[1].y) / 2);
                glEnd();
            } else {
                length = std::sqrt(xLen * xLen + yLen * yLen);
                intervalCount = (uint64_t)(length / interval);
                glBegin(GL_POINTS);
                    glColor3f((float) color.r / 255, (float) color.g / 255, (float) color.b / 255);
                    for (uint64_t i = 0; i < intervalCount; ++i) {
                        thr1 = {position[0].x + (position[1].x - position[0].x) / intervalCount / 2 * (i * 2 + 1), position[0].y + (position[1].y - position[0].y) / intervalCount / 2 * (i * 2 + 1)};
                        glVertex2f(thr1.x, thr1.y);
                    }
                glEnd();
            }
            break;
        case Type::LINE:
            glBegin(GL_LINES);
                glColor3f((float) color.r / 255, (float) color.g / 255, (float) color.b / 255);
                glVertex2f(position[0].x, position[0].y);
                glVertex2f(position[1].x, position[1].y);
            glEnd();
            break;
        case Type::LINE_WITH_POINT:
            glBegin(GL_LINES);
                glColor3f((float) color.r / 255, (float) color.g / 255, (float) color.b / 255);
                glVertex2f(position[0].x, position[0].y);
                glVertex2f(position[1].x, position[1].y);
            glEnd();
            if (interval == 0.f) {
                glBegin(GL_POINTS);
                    glColor3f((float) color.r / 255, (float) color.g / 255, (float) color.b / 255);
                    glVertex2f((position[0].x + position[1].x) / 2, (position[0].y + position[1].y) / 2);
                glEnd();
            } else {
                length = std::sqrt(xLen * xLen + yLen * yLen);
                intervalCount = (uint64_t)(length / interval);
                glBegin(GL_POINTS);
                    glColor3f((float) color.r / 255, (float) color.g / 255, (float) color.b / 255);
                    for (uint64_t i = 0; i < intervalCount; ++i) {
                        thr1 = {position[0].x + (position[1].x - position[0].x) / intervalCount / 2 * (i * 2 + 1), position[0].y + (position[1].y - position[0].y) / intervalCount / 2 * (i * 2 + 1)};
                        glVertex2f(thr1.x, thr1.y);
                    }
                glEnd();
            }
            break;
        case Type::DOTTED:
            if (interval == 0.f) {
                thr1 = {position[0].x + (position[1].x - position[0].x) / 3, position[0].y + (position[1].y - position[0].y) / 3};
                thr2 = {position[1].x - (position[1].x - position[0].x) / 3, position[1].y - (position[1].y - position[0].y) / 3};
                glBegin(GL_LINES);
                    glColor3f((float) color.r / 255, (float) color.g / 255, (float) color.b / 255);
                    glVertex2f(thr1.x, thr1.y);
                    glVertex2f(thr2.x, thr2.y);
                glEnd();
            } else {
                length = std::sqrt(xLen * xLen + yLen * yLen);
                intervalCount = (uint64_t)(length / interval);
                glBegin(GL_LINES);
                    glColor3f((float) color.r / 255, (float) color.g / 255, (float) color.b / 255);
                    for (uint64_t i = 0; i < intervalCount; ++i) {
                        thr1 = {position[0].x + (position[1].x - position[0].x) / intervalCount / 3 * (i * 3 + 1), position[0].y + (position[1].y - position[0].y) / intervalCount / 3 * (i * 3 + 1)};
                        thr2 = {position[0].x + (position[1].x - position[0].x) / intervalCount / 3 * (i * 3 + 2), position[0].y + (position[1].y - position[0].y) / intervalCount / 3 * (i * 3 + 2)};
                        glVertex2f(thr1.x, thr1.y);
                        glVertex2f(thr2.x, thr2.y);
                    }
                glEnd();
            }
            break;
        case Type::DOTTED_WITH_POINT:
            glBegin(GL_POINTS);
                glColor3f((float) color.r / 255, (float) color.g / 255, (float) color.b / 255);
                glVertex2f(position[0].x, position[0].y);
                glVertex2f(position[1].x, position[1].y);
            glEnd();
            thr1 = {position[0].x + (position[1].x - position[0].x) / 3, position[0].y + (position[1].y - position[0].y) / 3};
            thr2 = {position[1].x - (position[1].x - position[0].x) / 3, position[1].y - (position[1].y - position[0].y) / 3};
            glBegin(GL_LINES);
                glColor3f((float) color.r / 255, (float) color.g / 255, (float) color.b / 255);
                glVertex2f(thr1.x, thr1.y);
                glVertex2f(thr2.x, thr2.y);
            glEnd();
            break;
        default:
            glBegin(GL_LINES);
            break;
    }
    target.popGLStates();
    target.setActive(false);
}

Line::Line () {
    thickness = 1.f;
    interval = intervalMiss = 0.f;
    type = Type::LINE;
    color = sf::Color::Black;
    position[0] = position[1] = {0.f, 0.f};
}

Line::Line (const Line &line) {
    thickness = line.thickness;
    interval = line.interval;
    intervalMiss = line.intervalMiss;
    type = line.type;
    color = line.color;
    position[0] = line.position[0];
    position[1] = line.position[1];
}

Line::Line (Type type, sf::Color color, float thickness, float interval) : type(type), color(color), thickness(thickness), interval(interval) {
    intervalMiss = 0.f;
    position[0] = position[1] = {0.f, 0.f};
}

Line::~Line () {}

Line &Line::operator= (const Line &line) {
    thickness = line.thickness;
    interval = line.interval;
    intervalMiss = line.intervalMiss;
    type = line.type;
    color = line.color;
    position[0] = line.position[0];
    position[1] = line.position[1];
    return *this;
}