#include <GUI/Table.hpp>
#include <iostream>

const std::string Table::M_NAME = "Table";
const Table::GuiFabric<Table> Table::m_tableFabric(Table::M_NAME);

void Table::recalculateSize () {
    uint64_t charSize = m_text.getCharacterSize();
    m_VLengths.resize(m_content.size(), 0);
    m_HLengths.resize(m_content[0].size(), 0);
    for (uint64_t i = 0; i < m_content.size(); ++i) {
        for (uint64_t j = 0; j < m_content[i].size(); ++j) {
            m_text.setString(m_content[i][j]);
            m_HLengths[j] = std::max(m_text.getLocalBounds().width, m_HLengths[j]);
        }
    }
    for (uint64_t i = 0; i < m_content.size(); ++i) {
        for (uint64_t j = 0; j < m_content[i].size(); ++j) {
            m_text.setString(m_content[i][j]);
            float ratio = m_text.getLocalBounds().width / m_HLengths[j];
            uint64_t nLineCount = 1;
            if (ratio > 1.f) {
                nLineCount = (uint64_t)std::roundf(ratio + 0.5);
                for (uint64_t k = 0; k < nLineCount - 1; ++k) {
                    sf::String &str = m_content[i][j];
                    str.insert(str.getSize() / nLineCount + k, "\n");
                }
                //m_text.setString(m_content[i][j]);
            }
            m_VLengths[i] = std::max(m_VLengths[i], static_cast<float>(charSize * nLineCount));
        }
    }
    float maxHorizont = 0, maxVertical = 0;
    for (uint64_t i = 0; i < m_HLengths.size(); ++i) {
        maxHorizont += m_HLengths[i] + 10;
    }
    for (uint64_t i = 0; i < m_VLengths.size(); ++i) {
        maxVertical += m_VLengths[i] + 10;
    }
    setRenderSize(maxHorizont + 1, maxVertical + 1);
    m_size = {maxHorizont, maxVertical};
    m_needRedraw = true;
}

void Table::updateRender () {
    uint64_t charSize = m_text.getCharacterSize();
    m_render.clear(sf::Color::White);
    sf::RectangleShape rect;
    sf::Vertex lines[2];

    //rect.setFillColor(sf::Color(127, 127, 127, 127));
    rect.setFillColor(sf::Color::White);
    rect.setPosition(0, 0);
    lines[0].color = lines[1].color = sf::Color::Black;

    float offsetX = 0, offsetY = 0;
    for (uint64_t i = 0; i < m_content.size(); ++i) {
        if (m_grid && i % 2 == 0) {
            float prevPosY = rect.getPosition().y + rect.getSize().y;
            if (i > 0) {
                prevPosY += m_VLengths[i - 1] + 10;
            }
            rect.setSize({m_size.x, m_VLengths[i] + 10});
            rect.setPosition(0, prevPosY);
            m_render.draw(rect);
        }
        offsetX = 0;
        for (uint64_t j = 0; j < m_content[i].size(); ++j) {
            m_text.setString(m_content[i][j]);
            m_text.setPosition(offsetX + (j * 2 + 1) * 5, offsetY + (i * 2 + 1) * 5);
            m_render.draw(m_text);
            offsetX += m_HLengths[j];
        }
        offsetY += m_VLengths[i];
    }

    //cols and rows lines
    lines[0].position = {0, 0};
    lines[1].position = {m_size.x, 0};
    m_render.draw(lines, 2, sf::Lines);
    for (uint64_t i = 0; i < m_content.size(); ++i) {
        float prevPos = lines[0].position.y;
        lines[0].position = {0, m_VLengths[i] + prevPos + 10};
        lines[1].position = {m_size.x, m_VLengths[i] + prevPos + 10};
        m_render.draw(lines, 2, sf::Lines);
    }
    lines[0].position = {1, 0};
    lines[1].position = {1, m_size.y};
    m_render.draw(lines, 2, sf::Lines);
    for (uint64_t i = 0; i < m_content[0].size(); ++i) {
        float prevPos = lines[0].position.x;
        lines[0].position = {m_HLengths[i] + prevPos + 10, 0};
        lines[1].position = {m_HLengths[i] + prevPos + 10, m_size.y};
        m_render.draw(lines, 2, sf::Lines);
    }
    // if (m_field.getInteractionFocus()) {
    //     m_render.draw(m_field);
    //     m_render.draw(m_shape);
    // }
    m_render.display();
    std::cout << "BE Table\n";
    m_needRedraw = false;
}

void Table::updateContent () {
    if (m_field.getSignalStatus(Signal::VISUAL_CHANGE)) {
        handleSignal(Signal::VISUAL_CHANGE);
    }
    m_field.update();
    if (m_needRedraw || getSignalStatus(Signal::VISUAL_CHANGE)) {
        //recalculateSize();
        updateRender();
        m_needRedraw = false;
        handleSignal(Signal::VISUAL_CHANGE);
    }
    m_needUpdate = false;
}

void Table::interact (const sf::Event &event, const sf::Vector2f &mousePos) {
    static uint64_t cellX = 0, cellY = 0;
    if (event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == m_mouse) {
        if (contains(mousePos) && !m_field.getInteractionFocus()) {
            //m_needRedraw = true;
            //m_field.setInteractionFocus(true);
            cellX = 0, cellY = 0;
            float distX = 0, distY = 0;
            uint64_t charSize = m_field.getInputLine().getText().getCharacterSize();
            while (distX + m_HLengths[cellX] + 10 < mousePos.x || distY + m_VLengths[cellY] + 10 < mousePos.y) {
                if (distX + m_HLengths[cellX] + 10 < mousePos.x) {
                    distX += m_HLengths[cellX] + 10;
                    ++cellX;
                }
                if (distY + m_VLengths[cellY] + 10 < mousePos.y) {
                    distY += m_VLengths[cellY] + 10;
                    ++cellY;
                }
            }
            std::cout << "cell: " << cellX << " " << cellY << '\n';
            m_shape.setFillColor(sf::Color::Transparent);
            m_shape.setOutlineColor(sf::Color::Blue);
            m_shape.setOutlineThickness(-2.f);
            m_shape.setPosition(distX, distY);
            m_shape.setSize({m_HLengths[cellX] + 10, m_VLengths[cellY] + 10});

            m_field.getInputLine().setString(m_content[cellY][cellX]);
            m_field.setPosition(distX, distY);
            m_field.setSize({m_HLengths[cellX] + 10, m_VLengths[cellY] + 10});
            //m_field.getInputLine().setLimitOfSymbols((m_HLengths[cellX] + 10) / charSize);
            m_needUpdate = m_needRedraw = true;
            handleSignal(Signal::VISUAL_CHANGE);
            //m_field.setInteractionFocus(true);
        }
    }
    if (m_editable) {
        m_field.handleEvent(event, mousePos);
    }
    if (m_field.getSignalStatus(Signal::LOST_FOCUS)) {
        sf::String candidate = m_field.getInputLine().getString();
        if (candidate != m_content[cellY][cellX]) {
            m_content[cellY][cellX] = candidate;
            m_needRedraw = true;
        }
        handleSignal(Signal::VISUAL_CHANGE);
    }
    if (m_field.getSignalStatus(Signal::VISUAL_CHANGE)) {
        std::cout << "field change\n";
        handleSignal(Signal::VISUAL_CHANGE);
    }
}

void Table::drawVisible (sf::RenderTarget& target, sf::RenderStates states) const {
    target.draw(m_visibleContent, states);
    if (m_field.getInteractionFocus()) {
        target.draw(m_field, states);
        target.draw(m_shape, states);
    }
}

Table::Table () {
    m_grid = true;
    m_editable = false;
    m_field.setSignal(Signal::VISUAL_CHANGE, [this] () {
        this->handleSignal(Signal::VISUAL_CHANGE);
    });
}

Table::Table (const std::vector<std::vector<sf::String>> &content) : Table() {
    setContent(content);
}

Table::~Table () {}

void Table::setContent (const std::vector<std::vector<sf::String>> &content) {
    m_needUpdate = m_needRedraw = true;
    m_content = content;
    recalculateSize();
}

void Table::setCols (uint64_t cols) {
    m_needUpdate = m_needRedraw = true;
    for (uint64_t i = 0; i < m_content.size(); ++i) {
        m_content[i].resize(cols);
    }
}

void Table::setRows (uint64_t rows) {
    m_needUpdate = m_needRedraw = true;
    m_content.resize(rows);
}

void Table::addRow (uint64_t i, const std::vector<sf::String> &row) {
    m_needUpdate = m_needRedraw = true;
    m_content.insert(m_content.begin() + i, row);
    recalculateSize();
}

uint64_t Table::getRowCount () const {
    return m_content.size();
}

uint64_t Table::getColCount () const {
    if (!m_content.empty()) {
        return m_content[0].size();
    } else {
        return 0;
    }
}

void Table::setColSize (uint64_t i, float size) {
    m_needUpdate = m_needRedraw = true;
    m_HLengths[i] = size;
}

void Table::setRowSize (uint64_t i, float size) {
    m_needUpdate = m_needRedraw = true;
    m_VLengths[i] = size;
}

void Table::changeItem (uint64_t i, uint64_t j, const sf::String &item) {
    m_needUpdate = m_needRedraw = true;
    m_content[i][j] = item;
}

void Table::setGrid (bool grid) {
    m_needUpdate = m_needRedraw = true;
    m_grid = grid;
}

void Table::editStatus (bool status) {
    m_editable = status;
}

// void Table::update () {
//     Interactable::update();
//     if (m_needRedraw) {
//         updateRender();
//     }
// }

void Table::setSize (const sf::Vector2f &size) {
    m_needUpdate = m_needRedraw = true;
    float xDiff = (size.x - m_size.x) / m_HLengths.size(), yDiff = (size.y - m_size.y) / m_VLengths.size();
    for (uint64_t i = 0; i < m_HLengths.size(); ++i) {
        m_HLengths[i] += xDiff;
    }
    for (uint64_t i = 0; i < m_VLengths.size(); ++i) {
        m_VLengths[i] += yDiff;
    }
    m_size = size;
}

void Table::loadFromFile (std::ifstream &file) {
    uint64_t lenSize, cols, rows, strSize;
    std::basic_string<sf::Uint32> str;
    file.read(reinterpret_cast<char*>(&cols), sizeof(cols));
    if (cols > 0) {
        m_content.resize(cols);
    }
    for (uint64_t i = 0; i < cols; ++i) {
        file.read(reinterpret_cast<char*>(&rows), sizeof(rows));
        if (rows > 0) {
            m_content[i].resize(cols);
        }
        for (uint64_t j = 0; j < rows; ++j) {
            file.read(reinterpret_cast<char*>(&strSize), sizeof(strSize));
            if (strSize > 0) {
                str.resize(strSize);
                file.read(reinterpret_cast<char*>(&str[0]), sizeof(str[0]) * strSize);
                m_content[i][j] = str;
            }
        }
    }
    file.read(reinterpret_cast<char*>(&lenSize), sizeof(lenSize));
    if (lenSize > 0) {
        m_HLengths.resize(lenSize);
        file.read(reinterpret_cast<char*>(&m_HLengths[0]), sizeof(m_HLengths[0]) * lenSize);
    }
    file.read(reinterpret_cast<char*>(&lenSize), sizeof(lenSize));
    if (lenSize > 0) {
        m_VLengths.resize(lenSize);
        file.read(reinterpret_cast<char*>(&m_VLengths[0]), sizeof(m_VLengths[0]) * lenSize);
    }
    file.read(reinterpret_cast<char*>(&m_grid), sizeof(m_grid));
    file.read(reinterpret_cast<char*>(&m_editable), sizeof(m_editable));
    Interactable::loadFromFile(file);
    TextBased::loadFromFile(file);
    RenderTextureBased::loadFromFile(file);
}

void Table::saveToFile (std::ofstream &file) const {
    uint64_t lenSize, cols, rows, strSize;
    uint64_t size;
    cols = m_content.size();
    size = M_NAME.size();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    file.write(&M_NAME[0], sizeof(char) * size);
    file.write(reinterpret_cast<const char*>(&cols), sizeof(cols));
    for (uint64_t i = 0; i < cols; ++i) {
        rows = m_content[i].size();
        file.write(reinterpret_cast<const char*>(&rows), sizeof(rows));
        for (uint64_t j = 0; j < rows; ++j) {
            strSize = m_content[i][j].getSize();
            file.write(reinterpret_cast<const char*>(&strSize), sizeof(strSize));
            if (strSize > 0) {
                file.write(reinterpret_cast<const char*>(m_content[i][j].getData()), sizeof(m_content[i][j][0]) * strSize);
            }
        }
    }
    lenSize = m_HLengths.size();
    file.write(reinterpret_cast<const char*>(&lenSize), sizeof(lenSize));
    if (lenSize > 0) {
        file.write(reinterpret_cast<const char*>(&m_HLengths[0]), sizeof(m_HLengths[0]) * lenSize);
    }
    lenSize = m_VLengths.size();
    file.write(reinterpret_cast<const char*>(&lenSize), sizeof(lenSize));
    if (lenSize > 0) {
        file.write(reinterpret_cast<const char*>(&m_VLengths[0]), sizeof(m_VLengths[0]) * lenSize);
    }
    file.write(reinterpret_cast<const char*>(&m_grid), sizeof(m_grid));
    file.write(reinterpret_cast<const char*>(&m_editable), sizeof(m_editable));
    Interactable::saveToFile(file);
    TextBased::saveToFile(file);
    RenderTextureBased::saveToFile(file);
}

const std::vector<sf::String> &Table::operator[] (uint64_t i) const {
    return m_content[i];
}