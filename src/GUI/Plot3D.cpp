#include <GUI/Plot3D.hpp>
#include <iostream>

//recalculateSize
//setSize

// GraphInfo::GraphInfo () {}

// GraphInfo::GraphInfo (const std::vector<sf::Vector2f> &points, const sf::Color &color) : points(points), color(color) {}

// GraphInfo::~GraphInfo () {}

void Plot3D::createNormals () {
    // normals.resize(splitX * splitZ);
    // for (uint64_t i = 0; i < splitX; ++i) {
    //     for (uint64_t j = 0; j < splitZ; ++j) {
    //         normals[i * splitZ + j] = {0, 0, 0};
    //     }
    // }
    // // for (uint64_t i = 0; i < splitX - 1; ++i) {
    // //     for (uint64_t j = 0; j < splitZ - 1; ++j) {
    // //         const sf::Vector3f &p0 = plane[j], &p1 = plane[j + 1], &p2 = plane[(i + 1) * splitZ + j];
    // //         normals[i * splitZ + j] = {1, 1, 1};
    // //     }
    // // }

    // uint64_t k = 0;
    // for (uint64_t i = 0; i < splitX - 1; ++i) {
    //     for (uint64_t j = 0; j < splitZ - 1; ++j) {
    //         float a, b, c;
    //         const sf::Vector3f &p0 = plane[j], &p1 = plane[j + 1], &p2 = plane[(i + 1) * splitX + j];
    //         a = (p1.y - p0.y) * (p2.z - p0.z) - (p2.y - p0.y) * (p1.z - p0.z);
    //         b = (p1.x - p0.x) * (p2.z - p0.z) - (p2.x - p0.x) * (p1.z - p0.z);
    //         c = (p1.x - p0.x) * (p2.y - p0.y) - (p2.x - p0.x) * (p1.y - p0.y);
    //         normals[k] += {a, b, c};
    //         ++k;
    //     }
    //     normals[k] += normals[k - 1];
    //     ++k;
    // }
    // for (uint64_t i = 0; i < splitZ; ++i) {
    //     normals[k] += normals[(splitX - 2) * splitZ + i];
    //     ++k;
    // }
    //glGe
    normals.clear();
    for (uint64_t i = 0; i < splitX; ++i) {
        for (uint64_t j = 0; j < splitZ; ++j) {
            normals.push_back({1, 1, 1});
        }
    }
    uint64_t k = 0;
    for (uint64_t i = 0; i < splitX - 1; ++i) {
        for (uint64_t j = 0; j < splitZ - 1; ++j) {
            float a, b, c;
            const sf::Vector3f &p0 = plane[j], &p1 = plane[j + 1], &p2 = plane[(i + 1) * splitX + j];
            a = (p1.y - p0.y) * (p2.z - p0.z) - (p2.y - p0.y) * (p1.z - p0.z);
            b = (p1.x - p0.x) * (p2.z - p0.z) - (p2.x - p0.x) * (p1.z - p0.z);
            c = (p1.x - p0.x) * (p2.y - p0.y) - (p2.x - p0.x) * (p1.y - p0.y);
            normals[k] += {a, b, c};
            ++k;
        }
        normals[k] += normals[k - 1];
        ++k;
    }
    for (uint64_t i = 0; i < splitZ; ++i) {
        normals[k] += normals[(splitX - 2) * splitZ + i];
        ++k;
    }
}

void Plot3D::recalculateSize () {

}

void Plot3D::updateRender () {
    uint64_t k = 0;
    m_render.clear(m_clearColor);
    m_render.setActive(true);
    // Configure the viewport (the same size as the window)
    glViewport(0, 0, m_render.getSize().x, m_render.getSize().y);
    
    //m_render.pushGLStates();
    //glBindBuffer( GL_ARRAY_BUFFER, 0 );
    //glDisableVertexAttribArray( 0 );
    // m_render.resetGLStates();
    // for (uint64_t i = 0; i < 5; ++i) {
        
    //     m_text.setString(std::to_string(i));
    //     m_text.setPosition(startX * scaleX + 10 * (i * 5) + 15, startY * scaleY + 10 * i);
    //     //m_render.pushGLStates();
    //     m_render.draw(m_text);
    //     //m_render.popGLStates();

    // }
    //m_render.popGLStates();

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LINE_SMOOTH);
    glDepthMask(GL_TRUE);
    glDepthFunc(GL_LESS);
    //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    //glClearDepth(1.f);
    // Setup a perspective projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    double ratio = static_cast<double>(m_render.getSize().x) / m_render.getSize().y;
    glFrustum(-ratio, ratio, -1.0, 1.0, 1.0, 500.0);
    glTranslatef(0, 0, -100.f);
    // Bind the texture
    //glEnable(GL_TEXTURE_2D);
    //glMatrixMode(GL_MODELVIEW);
    //glLoadIdentity();
    //glEnable(GL_CULL_FACE);
    //glCullFace(GL_FRONT);
    glEnable(GL_LIGHT0);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_NORMALIZE);
    // glEnable(GL_DEPTH_TEST);
    // glDepthMask(GL_TRUE);
    // glClearDepth(1.f);
    //glTranslatef(0, 0, -100.f);
    glPushMatrix();

    float lightPos[] = {0, 0, 100, 0};
    glLightfv(GL_LIGHT0, GL_POSITION, lightPos);
    glRotatef(rotateX, 1.0f, 0.0f, 0.0f);
    glRotatef(rotateY, 0.0f, 1.0f, 0.0f);
    glScalef(m_scale, m_scale, m_scale);
    glLineWidth(2.f);
    //отрисовка поверхности
    if (m_planeColor != sf::Color::Transparent) {
        if (m_needLight) {
            glEnable(GL_LIGHTING);
        }
        glBegin(GL_QUADS);
            glColor3f((float) m_planeColor.r / 255, (float) m_planeColor.g / 255, (float) m_planeColor.b / 255);
            for (uint64_t i = 0; i < splitX - 1; ++i) {
                for (uint64_t j = 0; j < splitZ - 1; ++j) {
                    k = i * splitZ + j;
                    glNormal3f(normals[k].x, normals[k].y, normals[k].z);
                    glVertex3f(plane[k].x * scaleX, plane[k].y * scaleY, plane[k].z * scaleZ);
                    k = i * splitZ + j + 1;
                    glNormal3f(normals[k].x, normals[k].y, normals[k].z);
                    glVertex3f(plane[k].x * scaleX, plane[k].y * scaleY, plane[k].z * scaleZ);
                    k = (i + 1) * splitZ + j + 1;
                    glNormal3f(normals[k].x, normals[k].y, normals[k].z);
                    glVertex3f(plane[k].x * scaleX, plane[k].y * scaleY, plane[k].z * scaleZ);
                    k = (i + 1) * splitZ + j;
                    glNormal3f(normals[k].x, normals[k].y, normals[k].z);
                    glVertex3f(plane[k].x * scaleX, plane[k].y * scaleY, plane[k].z * scaleZ);
                }
            }
        glEnd();
        if (m_needLight) {
            glDisable(GL_LIGHTING);
        }
    }
    glDisable(GL_DEPTH_TEST);

    //отрисовка сетки поверхности
    if (m_lineColor != sf::Color::Transparent) {
        for (uint64_t i = 0; i < splitX; ++i) {
            glBegin(GL_LINE_STRIP);
                glColor3f((float) m_lineColor.r / 255, (float) m_lineColor.g / 255, (float) m_lineColor.b / 255);
                for (uint64_t j = 0; j < splitZ; ++j) {
                    k = i * splitZ + j;
                    glVertex3f(plane[k].x * scaleX, plane[k].y * scaleY, plane[k].z * scaleZ);
                }
            glEnd();
        }
        for (uint64_t i = 0; i < splitX; ++i) {
            glBegin(GL_LINE_STRIP);
                glColor3f((float) m_lineColor.r / 255, (float) m_lineColor.g / 255, (float) m_lineColor.b / 255);
                for (uint64_t j = 0; j < splitZ; ++j) {
                    k = j * splitX + i;
                    glVertex3f(plane[k].x * scaleX, plane[k].y * scaleY, plane[k].z * scaleZ);
                }
            glEnd();
        }
    }
    sf::Color axisColor = sf::Color(127, 127, 127);
    float lilStep, bigStep;

    //отрисовка координат
    glBegin(GL_LINES);
        glColor3f((float) axisColor.r / 255, (float) axisColor.g / 255, (float) axisColor.b / 255);

        glVertex3f(startX * scaleX, startY * scaleY, startZ * scaleZ);
        glVertex3f(endX * scaleX, startY * scaleY, startZ * scaleZ);
        lilStep = (endX - startX) * scaleX / 25;
        bigStep = lilStep * 5;
        for (uint64_t i = 0; i < 5; ++i) {
            for (uint64_t j = 0; j < 5; ++j) {
                glVertex3f(startX * scaleX + lilStep * (i * 5 + j), startY * scaleY, startZ * scaleZ);
                glVertex3f(startX * scaleX + lilStep * (i * 5 + j), startY * scaleY + 5, startZ * scaleZ);
                glVertex3f(startX * scaleX + lilStep * (i * 5 + j), startY * scaleY, startZ * scaleZ);
                glVertex3f(startX * scaleX + lilStep * (i * 5 + j), startY * scaleY, startZ * scaleZ + 5);
            }
            glVertex3f(startX * scaleX + lilStep * ((i + 1) * 5), startY * scaleY, startZ * scaleZ);
            glVertex3f(startX * scaleX + lilStep * ((i + 1) * 5), startY * scaleY + 10, startZ * scaleZ);
            glVertex3f(startX * scaleX + lilStep * ((i + 1) * 5), startY * scaleY, startZ * scaleZ);
            glVertex3f(startX * scaleX + lilStep * ((i + 1) * 5), startY * scaleY, startZ * scaleZ + 10);

            glVertex3f(startX * scaleX + lilStep * ((i + 1) * 5), startY * scaleY, startZ * scaleZ);
            glVertex3f(startX * scaleX + lilStep * ((i + 1) * 5), startY * scaleY, endZ * scaleZ);
            //m_text.setString("1");
            //m_text.setPosition(startX * scaleX + lilStep * (i * 5), startY * scaleY);
            //m_render.pushGLStates();
            //m_render.draw(m_text);
            //m_render.popGLStates();

        }

        glVertex3f(startX * scaleX, startY * scaleY, startZ * scaleZ);
        glVertex3f(startX * scaleX, endY * scaleY, startZ * scaleZ);
        lilStep = (endY - startY) * scaleY / 25;
        bigStep = lilStep * 5;
        for (uint64_t i = 0; i < 5; ++i) {
            for (uint64_t j = 0; j < 5; ++j) {
                glVertex3f(startX * scaleX, startY * scaleY + lilStep * (i * 5 + j), startZ * scaleZ);
                glVertex3f(startX * scaleX + 5, startY * scaleY + lilStep * (i * 5 + j), startZ * scaleZ);
                glVertex3f(startX * scaleX, startY * scaleY + lilStep * (i * 5 + j), startZ * scaleZ);
                glVertex3f(startX * scaleX, startY * scaleY + lilStep * (i * 5 + j), startZ * scaleZ + 5);
            }
            glVertex3f(startX * scaleX, startY * scaleY + lilStep * ((i + 1) * 5), startZ * scaleZ);
            glVertex3f(startX * scaleX + 10, startY * scaleY + lilStep * ((i + 1) * 5), startZ * scaleZ);
            glVertex3f(startX * scaleX, startY * scaleY + lilStep * ((i + 1) * 5), startZ * scaleZ);
            glVertex3f(startX * scaleX, startY * scaleY + lilStep * ((i + 1) * 5), startZ * scaleZ + 10);
        }

        glVertex3f(startX * scaleX, startY * scaleY, startZ * scaleZ);
        glVertex3f(startX * scaleX, startY * scaleY, endZ * scaleZ);
        lilStep = (endZ - startZ) * scaleZ / 25;
        bigStep = lilStep * 5;
        for (uint64_t i = 0; i < 5; ++i) {
            for (uint64_t j = 0; j < 5; ++j) {
                glVertex3f(startX * scaleX, startY * scaleY, startZ * scaleZ + lilStep * (i * 5 + j));
                glVertex3f(startX * scaleX + 5, startY * scaleY, startZ * scaleZ + lilStep * (i * 5 + j));
                glVertex3f(startX * scaleX, startY * scaleY, startZ * scaleZ + lilStep * (i * 5 + j));
                glVertex3f(startX * scaleX, startY * scaleY + 5, startZ * scaleZ + lilStep * (i * 5 + j));
            }
            glVertex3f(startX * scaleX, startY * scaleY, startZ * scaleZ + lilStep * ((i + 1) * 5));
            glVertex3f(startX * scaleX + 10, startY * scaleY, startZ * scaleZ + lilStep * ((i + 1) * 5));
            glVertex3f(startX * scaleX, startY * scaleY, startZ * scaleZ + lilStep * ((i + 1) * 5));
            glVertex3f(startX * scaleX, startY * scaleY + 10, startZ * scaleZ + lilStep * ((i + 1) * 5));

            glVertex3f(startX * scaleX, startY * scaleY, startZ * scaleZ + lilStep * ((i + 1) * 5));
            glVertex3f(endX * scaleX, startY * scaleY, startZ * scaleZ + lilStep * ((i + 1) * 5));
        }
    glEnd();
    //glPopMatrix();
    //m_render.pushGLStates();
    
    //m_render.popGLStates();
    m_render.display();
    m_render.setActive(false);
    m_needRedraw = false;
}

void Plot3D::updateContent () {
    if (m_needRedraw) {
        updateRender();
        m_needRedraw = false;
    }
}

void Plot3D::interact (const sf::Event &event, const sf::Vector2f &mousePos) {
    auto newMousePos = mousePos;
    if (contains(mousePos)) {
        if (newMousePos.x > 0 && newMousePos.x < 400 && newMousePos.y > 0 && newMousePos.y < 400) {
            currPos.x = newMousePos.x / scaleX + startX;
            currPos.y = -newMousePos.y / scaleY + endY;
            m_text.setString("(" + toString(currPos.x, 2) + "; " + toString(currPos.y, 2) + ")");
        }
        if (event.type == sf::Event::MouseWheelMoved) {
            m_scale += 0.01 * event.mouseWheel.delta;
            m_needRedraw = true;
        }
    }
    if (m_interactFocus) {
        if (event.type == sf::Event::MouseMoved && sf::Mouse::isButtonPressed(sf::Mouse::Left)) {
            rotateX += (newMousePos.y - oldMousePos.y) * 0.2f;
            rotateY += (newMousePos.x - oldMousePos.x) * 0.2f;
            //rotateX = newMousePos.y * 0.5f;
            //rotateY = newMousePos.x * 0.5f;
            // float distX = endX - startX, distY = endY - startY;
            // startX += (oldMousePos.x - newMousePos.x) / scaleX;
            // endX += (oldMousePos.x - newMousePos.x) / scaleX;
            // startY += (newMousePos.y - oldMousePos.y) / scaleY;
            // endY += (newMousePos.y - oldMousePos.y) / scaleY;
            m_needRedraw = true;
        }
        if (event.type == sf::Event::MouseWheelMoved) {
            m_scale += 0.01 * event.mouseWheel.delta;
            m_needRedraw = true;
        }
        if (event.type == sf::Event::KeyPressed) {
            if (event.key.code == sf::Keyboard::Key::Left) {

            }
            m_needRedraw = true;
        }
        oldMousePos = mousePos;
    }
}

void Plot3D::drawVisible (sf::RenderTarget& target, sf::RenderStates states) const {
    //target.draw(border, states);
    target.draw(m_visibleContent, states);
    // if (showPoint) {
    //     target.draw(m_text, states);
    // }
}

sf::Color Plot3D::getRandomColor () {
    uint32_t num = static_cast<uint32_t>(std::rand());
    //num = num & (((255 << 8 + 255) << 8 + 255) << 8);
    num = (num << 8) + 255;
    return sf::Color(num);
}

void Plot3D::setAxisBorders () {
    for (uint64_t i = 0; i < plane.size(); ++i) {
        startX = std::min(startX, plane[i].x);
        endX = std::max(plane[i].x, endX);
        startY = std::min(startY, plane[i].y);
        endY = std::max(plane[i].y, endY);
        startZ = std::min(startZ, plane[i].z);
        endZ = std::max(plane[i].z, endZ);
    }
}

void Plot3D::setAxisScales () {
    scaleX = 400.f / (endX - startX);
    scaleY = 400.f / (endY - startY);
    scaleZ = 400.f / (endZ - startZ);
}

// void Plot3D::addOffset (const sf::Vector2f &off) {
//     lines[0].position += off;
//     lines[1].position += off;
// }

Plot3D::Plot3D () {
    m_size = {400, 400};
    setRenderSize(400, 400);
    rotateX = rotateY = 0;
    m_scale = 1.0;
    //m_render.create(800, 500);
    //visibleContent.setTexture(m_render.getTexture());
    // legendShape.setFillColor(sf::Color::White);
    // legendShape.setOutlineColor(sf::Color::Black);
    // legendShape.setOutlineThickness(1.f);
    m_needTriggerPressed = true;
    m_needTriggerInside = false;
    m_pressOutsideToLoseFocus = true;
    splitX = splitY = splitZ = 5;
    m_needGrid = m_needLegend = m_showPoints = m_needLight = true;
}

Plot3D::Plot3D (const std::function<float (float, float)> &func, uint64_t iterLimit, float x1, float x2, float z1, float z2) : Plot3D() {
    setPlane(func, iterLimit, x1, x2, z1, z2);
    m_needUpdate = m_needRedraw = true;
}

Plot3D::Plot3D (const std::vector<sf::Vector3f> &points, uint64_t iterX, uint64_t iterZ) : Plot3D() {
    setPlane(points, iterX, iterZ);
    m_needUpdate = m_needRedraw = true;
}

Plot3D::~Plot3D () {}

void Plot3D::setPlane (const std::function<float (float, float)> &func, uint64_t iterLimit, float x1, float x2, float z1, float z2) {
    float xStep = (x2 - x1) / iterLimit, zStep = (z2 - z1) / iterLimit;
    plane.resize(iterLimit * iterLimit);
    splitX = iterLimit;
    splitZ = iterLimit;
    for (uint64_t i = 0; i < iterLimit; ++i) {
        for (uint64_t j = 0; j < iterLimit; ++j) {
            plane[i * iterLimit + j] = {xStep * i, func(xStep * i, zStep * j), zStep * j};
        }
    }
    createNormals();
    setAxisBorders();
    setAxisScales();
}

void Plot3D::setPlane (const std::vector<sf::Vector3f> &points, uint64_t iterX, uint64_t iterZ) {
    plane = points;
    splitX = iterX;
    splitZ = iterZ;
    createNormals();
    setAxisBorders();
    setAxisScales();
}

void Plot3D::setXBorders (float x1, float x2) {
    startX = x1;
    endX = x2;
}

void Plot3D::setYBorders (float y1, float y2) {
    startY = y1;
    endY = y2;
}

void Plot3D::setZBorders (float z1, float z2) {
    startZ = z1;
    endZ = z2;
}

// void Plot3D::addLegend (uint64_t i, const sf::String &name) {
//     legend[i] = name;
//     updateContent();
// }

// void Plot3D::addLegends (const std::vector<sf::String> &names) {
//     legend = names;
//     updateContent();
// }

void Plot3D::setLineColor (const sf::Color &color) {
    m_lineColor = color;
    m_needUpdate = true;
}

void Plot3D::setPlaneColor (const sf::Color &color) {
    m_planeColor = color;
    m_needUpdate = true;
}

void Plot3D::clear () {
    startX = startY = INFINITY;
    endX = endY = -INFINITY;
    startZ = endZ = INFINITY;
    plane.clear();
    normals.clear();
    m_needUpdate = m_needRedraw = true;
}

void Plot3D::setSize (const sf::Vector2f &v) {}