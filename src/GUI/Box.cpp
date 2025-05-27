#include <GUI/Box.hpp>
#include <iostream>

const std::string Box::M_NAME = "Box";
//Box::BoxStaticConstructor Box::m_boxStaticConsturctor(Box::M_NAME);
const Box::GuiFabric<Box> Box::m_boxFabric(Box::M_NAME);

void Box::recalculateSize () {
    if (m_content.empty()) {
        return;
    }
    float width = 0, height = 0;
    m_maxWidth = m_maxHeight = 0;
    for (const auto &el : m_content) {
        // if (el.second->getType() == GUIType::CONTAINER) {
        //     auto container = std::static_pointer_cast<Box>(el.second);
        //     container->recalculateSize();
        //     container->realign();
        // }
        auto size = el.second->getSize();
        width += size.x;
        height += size.y;
        m_maxWidth = std::max(m_maxWidth, size.x);
        m_maxHeight = std::max(m_maxHeight, size.y);
    }
    if (m_orientation == Orientation::HORIZONTAL) {
        width += (m_content.size() - 1) * m_offset;
        height = m_maxHeight;
    } else {
        width = m_maxWidth;
        height += (m_content.size() - 1) * m_offset;
    }
    m_size = {width, height};
    m_needRealign = true;
    //std::cout << "size: " << m_size.x << " " << m_size.y << "\n";
}

void Box::realign () {
    if (m_content.size() <= 1) {
        return;
    }
    //std::cout << "realigning\n";
    using MapIterator = decltype(m_content)::iterator;
    if (m_needUniteSize) {
        if (m_orientation == Orientation::HORIZONTAL) {
            for (auto &el : m_content) {
                el.second->setSize(el.second->getSize().x, m_size.y);
                //std::cout << "new size of " << el.first << ": " << el.second->getSize().x << " " << el.second->getSize().y << "\n";
            }
        } else {
            for (auto &el : m_content) {
                el.second->setSize(m_size.x, el.second->getSize().y);
                std::cout << "new size of " << el.first << ": " << el.second->getSize().x << " " << el.second->getSize().y << "\n";
            }
        }
    }
    std::function<float (std::shared_ptr<GUIElement> &el, uint64_t i)> xStep, yStep;
    if (m_orientation == Orientation::HORIZONTAL) {
        xStep = [&] (std::shared_ptr<GUIElement> &el, uint64_t i) -> float {
            if (i == 0) {
                return 0;
            }
            auto &curr = m_content[m_names[i - 1]];
            return curr->getPosition().x + curr->getSize().x + m_offset;
        };
        switch (m_alignment) {
            case Alignment::UP:
                yStep = [&] (std::shared_ptr<GUIElement> &el, uint64_t i) {
                    return 0;
                };
                break;
            case Alignment::DOWN:
                yStep = [&] (std::shared_ptr<GUIElement> &el, uint64_t i) {
                    return m_size.y - el->getSize().y;
                };
                break;
            case Alignment::CENTER:
                yStep = [&] (std::shared_ptr<GUIElement> &el, uint64_t i) {
                    return (m_size.y - el->getSize().y) / 2;
                };
                break;
            default:
                yStep = [&] (std::shared_ptr<GUIElement> &el, uint64_t i) {
                    return 0;
                };
                break;
        }
    } else {
        yStep = [&] (std::shared_ptr<GUIElement> &el, uint64_t i) -> float {
            if (i == 0) {
                return 0;
            }
            auto &curr = m_content[m_names[i - 1]];
            return curr->getPosition().y + curr->getSize().y + m_offset;
        };
        switch (m_alignment) {
            case Alignment::LEFT:
                xStep = [&] (std::shared_ptr<GUIElement> &el, double offset) {
                    return 0;
                };
                break;
            case Alignment::RIGHT:
                xStep = [&] (std::shared_ptr<GUIElement> &el, double offset) {
                    return m_size.x - el->getSize().x;
                };
                break;
            case Alignment::CENTER:
                xStep = [&] (std::shared_ptr<GUIElement> &el, double offset) {
                    return (m_size.x - el->getSize().x) / 2;
                };
                break;
            default:
                xStep = [&] (std::shared_ptr<GUIElement> &el, double offset) {
                    return 0;
                };
                break;
        }
    }
    for (uint64_t i = 0; i < m_names.size(); ++i) {
        auto &curr = m_content[m_names[i]];
        curr->setPosition(xStep(curr, i), yStep(curr, i));
    }
    // while (it != itEnd) {
    //     std::cout << "q";
    //     it->second->setPosition(xStep(it, itBegin), yStep(it, itBegin));
    //     //it = std::next(it);
    //     ++it;
    // }
    //itEnd->second->setPosition(xStep(itEnd, itBegin), yStep(itEnd, itBegin));
    handleSignal(Signal::VISUAL_CHANGE);
    m_needRealign = false;
}

void Box::updateContent () {
    GUIContainer::updateContent();
    if (m_needRealign) {
        //recalculateSize();
        realign();
    }

}

void Box::drawVisible (sf::RenderTarget& target, sf::RenderStates states) const {
    GUIContainer::drawVisible(target, states);
}

Box::Box () {
    m_type = GUIType::CONTAINER;
    m_maxWidth = m_maxHeight = 0;
    m_offset = 20;
    m_alignment = Alignment::CENTER;
    m_orientation = Orientation::VERTICAL;
    //m_border.setFillColor(sf::Color::Transparent);
    //m_border.setOutlineColor(sf::Color::Black);
    //m_border.setOutlineThickness(-2.f);
    //m_needBorder = true;
    m_needRealign = true;
    m_needUniteSize = true;
    m_saveOffset = false;
    m_size = {100, 50};
}

Box::Box (uint64_t offset, Alignment alignment, Orientation orientation) : Box() {
    m_offset = offset;
    m_alignment = alignment;
    m_orientation = orientation;
}

Box::~Box () {}

void Box::setOffset (uint64_t offset) {
    m_needRealign = true;
    m_offset = offset;
}

void Box::setAlignment (Alignment alignment) {
    m_needRealign = true;
    m_alignment = alignment;
}

void Box::setOrientation (Orientation orientation) {
    m_needRealign = true;
    m_orientation = orientation;
}

void Box::needUniteSize (bool status) {
    m_needRealign = true;
    m_needUniteSize = status;
}

void Box::saveOffset (bool status) {
    m_needRealign = true;
    m_saveOffset = status;
}

void Box::setSize (const sf::Vector2f &size) {
    if (m_content.empty()) {
        m_size = size;
        return;
    }
    if (m_content.size() == 1) {
        // if (m_needUniteSize) {
        //     for (auto &el : m_content) {
        //         el.second->setSize(size);
        //     }
        // }
        m_size = size;
        if (m_needUniteSize) {
            if (m_orientation == Orientation::HORIZONTAL) {
                for (auto &el : m_content) {
                    el.second->setSize(el.second->getSize().x, m_size.y);
                }
            } else {
                for (auto &el : m_content) {
                    el.second->setSize(m_size.x, el.second->getSize().y);
                }
            }
        }
        //m_needRealign = true;
        return;
    }
    uint64_t count = m_content.size() - 1;
    float freeSpaceX = size.x - m_size.x + m_offset * count;
    float freeSpaceY = size.y - m_size.y + m_offset * count;
    float newOffset = 0;
    // if (m_needUniteSize) {
    //     count = 0;
    //     for (auto &el : m_content) {
    //         if (el.second->isResizeable()) {
    //             ++count;
    //         }
    //     }
    //     //freeSpaceX = freeSpaceX / (count + 1);//(freeSpaceX - m_offset * (m_content.size() - 1)) / count;
    //     //freeSpaceY = freeSpaceY / (count + 1);//(freeSpaceY - m_offset * (m_content.size() - 1)) / count;
    //     //freeSpaceX = size.x - m_size.x + m_offset * count;
    //     //freeSpaceY = size.y - m_size.y + m_offset * count;
    //     if (m_orientation == Orientation::HORIZONTAL) {
    //         float yAppend = freeSpaceY / (count + 1);
    //         for (auto &el : m_content) {
    //             if (el.second->isResizeable()) {
    //                 auto prevSize = el.second->getSize();
    //                 el.second->setSize(prevSize + sf::Vector2f(0, yAppend));
    //             }
    //         }
    //     } else {
    //         float xAppend = freeSpaceX / (count + 1);
    //         for (auto &el : m_content) {
    //             if (el.second->isResizeable()) {
    //                 auto prevSize = el.second->getSize();
    //                 el.second->setSize(prevSize + sf::Vector2f(xAppend, 0));
    //             }
    //         }
    //     }
    // }
    newOffset = m_orientation == Orientation::HORIZONTAL ? freeSpaceX / count : freeSpaceY / count;
    if (newOffset > M_MINIMAL_OFFSET) {
        m_offset = newOffset;
    } else {
        std::cout << "Box: bad size\n";
    }
    m_size = size;
    m_needRealign = true;
}

void Box::loadFromFile (std::ifstream &file) {
    file.read(reinterpret_cast<char*>(&m_maxWidth), sizeof(m_maxWidth));
    file.read(reinterpret_cast<char*>(&m_maxHeight), sizeof(m_maxHeight));
    file.read(reinterpret_cast<char*>(&m_offset), sizeof(m_offset));
    file.read(reinterpret_cast<char*>(&m_alignment), sizeof(m_alignment));
    file.read(reinterpret_cast<char*>(&m_orientation), sizeof(m_orientation));
    file.read(reinterpret_cast<char*>(&m_needRealign), sizeof(m_needRealign));
    file.read(reinterpret_cast<char*>(&m_needUniteSize), sizeof(m_needUniteSize));
    file.read(reinterpret_cast<char*>(&m_saveOffset), sizeof(m_saveOffset));
    GUIContainer::loadFromFile(file);
}

void Box::saveToFile (std::ofstream &file) const {
    uint64_t size = M_NAME.size();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    file.write(&M_NAME[0], sizeof(char) * size);
    file.write(reinterpret_cast<const char*>(&m_maxWidth), sizeof(m_maxWidth));
    file.write(reinterpret_cast<const char*>(&m_maxHeight), sizeof(m_maxHeight));
    file.write(reinterpret_cast<const char*>(&m_offset), sizeof(m_offset));
    file.write(reinterpret_cast<const char*>(&m_alignment), sizeof(m_alignment));
    file.write(reinterpret_cast<const char*>(&m_orientation), sizeof(m_orientation));
    file.write(reinterpret_cast<const char*>(&m_needRealign), sizeof(m_needRealign));
    file.write(reinterpret_cast<const char*>(&m_needUniteSize), sizeof(m_needUniteSize));
    file.write(reinterpret_cast<const char*>(&m_saveOffset), sizeof(m_saveOffset));
    GUIContainer::saveToFile(file);
}