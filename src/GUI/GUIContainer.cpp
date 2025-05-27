#include <GUI/GUIContainer.hpp>

void GUIContainer::updateContent () {
    for (auto &el : m_content) {
        if (el.second->getType() == GUIType::INTERACTABLE || el.second->getType() == GUIType::CONTAINER) {
            auto interactable = std::static_pointer_cast<Interactable>(el.second);
            interactable->update();
        }
    }
    if (getSignalStatus(Signal::ADDED_ITEM)) {
        std::cout << "GuiContainer: recalculating size\n";
        recalculateSize();
    }
}

void GUIContainer::interact (const sf::Event &event, const sf::Vector2f &mousePos) {
    m_interactFocus = false;
    for (auto &el : m_content) {
        if (el.second->getType() == GUIType::INTERACTABLE || el.second->getType() == GUIType::CONTAINER) {
            auto interactable = std::static_pointer_cast<Interactable>(el.second);
            interactable->handleEvent(event, mousePos);
            // if (interactable->getSignalStatus(Signal::VISUAL_CHANGE)) {
            //     std::cout << "b";
            //     handleSignal(Signal::VISUAL_CHANGE);
            //     //m_needRecalculateSize = true;
            // }
            if (interactable->getInteractionFocus()) {
                auto it = std::find(m_names.begin(), m_names.end(), el.first);
                m_activeObject = std::distance(m_names.begin(), it);
                m_interactFocus = true;
            }
        }
    }
}

void GUIContainer::drawVisible (sf::RenderTarget& target, sf::RenderStates states) const {
    if (!m_content.empty()) {
        for (uint64_t i = 0; i < m_names.size(); ++i) {
            if (i != m_activeObject) {
                target.draw(*m_content.at(m_names[i]), states);
            }
        }
        if (m_activeObject < m_names.size()) {
            target.draw(*m_content.at(m_names[m_activeObject]), states);
        }
    }
}

GUIContainer::GUIContainer () {
    m_activeObject = -1;
    //m_needRecalculateSize = false;
}

GUIContainer::~GUIContainer () {}

void GUIContainer::add (const std::string &name, const std::string &gui_el) {
    auto &fabric = getFabric();
    std::shared_ptr<GUIElement> element;
    if (fabric.count(gui_el) != 1) {
        throw std::runtime_error("GuiContainer::add: unknown element \"" + gui_el + "\"");
    }
    element = fabric[gui_el]();
    add(name, std::move(element));
}

void GUIContainer::add (const std::string &name, std::shared_ptr<GUIElement> &&element) {
    if (name.empty()) {
        throw std::runtime_error("GuiContainer::add: empty key not allowed");
    }
    if (m_content.count(name) != 0) {
        throw std::runtime_error("GuiContainer::add: item with key \"" + name + "\" already presented");
    }
    if (element->getType() == GUIType::INTERACTABLE || element->getType() == GUIType::CONTAINER) {
        auto interactable = std::static_pointer_cast<Interactable>(element);
        interactable->setSignal(Signal::VISUAL_CHANGE, [this] () {
            this->handleSignal(Signal::VISUAL_CHANGE);
        });
    }
    //m_content[name] = std::move(element);
    m_content[name] = element;
    m_names.push_back(name);
    //m_needRecalculateSize = true;
    handleSignal(Signal::VISUAL_CHANGE);
    handleSignal(Signal::ADDED_ITEM);
}

void GUIContainer::remove (const std::string &name) {
    m_content.erase(name);
    auto it = std::find(m_names.begin(), m_names.end(), name);
    if (it == m_names.begin() + m_activeObject) {
        m_names.erase(it);
        m_activeObject = -1;
    }
    //m_needRecalculateSize = true;
    handleSignal(Signal::VISUAL_CHANGE);
    handleSignal(Signal::ADDED_ITEM);
}

// template <class T>
// std::shared_ptr<T> &GUIContainer::get (const std::string &name) {
//     return std::dynamic_pointer_cast<T>(m_content.at(name));
// }

std::shared_ptr<GUIElement> &GUIContainer::getElement (const std::string &name) {
    return m_content[name];
}

const std::shared_ptr<GUIElement> &GUIContainer::getElement (const std::string &name) const {
    return m_content.at(name);
}

std::shared_ptr<GUIContainer> GUIContainer::getContainer (const std::string &name) {
    return std::dynamic_pointer_cast<GUIContainer>(m_content.at(name));
}

std::shared_ptr<const GUIContainer> GUIContainer::getContainer (const std::string &name) const {
    return std::dynamic_pointer_cast<const GUIContainer>(m_content.at(name));
}

std::shared_ptr<Interactable> GUIContainer::getInteractable (const std::string &name) {
    return std::dynamic_pointer_cast<Interactable>(m_content.at(name));
}

std::shared_ptr<const Interactable> GUIContainer::getInteractable (const std::string &name) const {
    return std::dynamic_pointer_cast<const Interactable>(m_content.at(name));
}

uint64_t GUIContainer::getContentSize () const {
    return m_content.size();
}

void GUIContainer::forceRecalculateSize () {
    //m_needRecalculateSize = true;
    handleSignal(Signal::VISUAL_CHANGE);
}

std::shared_ptr<GUIContainer> GUIContainer::at (const std::string &name) {
    return std::dynamic_pointer_cast<GUIContainer>(m_content[name]);
}

std::shared_ptr<const GUIContainer> GUIContainer::at (const std::string &name) const {
    return std::dynamic_pointer_cast<const GUIContainer>(m_content.at(name));
}

std::shared_ptr<GUIContainer> GUIContainer::operator[] (const std::string &name) {
    return std::dynamic_pointer_cast<GUIContainer>(m_content[name]);
}

std::shared_ptr<const GUIContainer> GUIContainer::operator[] (const std::string &name) const {
    return std::dynamic_pointer_cast<const GUIContainer>(m_content.at(name));
}

// std::shared_ptr<GUIContainer> GUIContainer::operator[] (uint64_t i) {
//     return std::dynamic_pointer_cast<GUIContainer>(m_content[m_names[i]]);
// }

// const std::shared_ptr<GUIContainer> GUIContainer::operator[] (uint64_t i) const {
//     return std::dynamic_pointer_cast<GUIContainer>(m_content.at(m_names[i]));
// }

void GUIContainer::loadFromFile (std::ifstream &file) {
    uint64_t strSize, contentSize;
    std::string id, name;
    file.read(reinterpret_cast<char*>(&contentSize), sizeof(contentSize));
    for (uint64_t i = 0; i < contentSize; ++i) {
        file.read(reinterpret_cast<char*>(&strSize), sizeof(strSize));
        id.resize(strSize);
        file.read(&id[0], sizeof(char) * strSize);
        file.read(reinterpret_cast<char*>(&strSize), sizeof(strSize));
        name.resize(strSize);
        file.read(&name[0], sizeof(char) * strSize);
        std::shared_ptr<GUIElement> element;
        //if (name != "Button") {
        //    element = GuiFabric::createFromString(name);
        //} else {
        auto &fabric = getFabric();
        if (fabric.count(name) != 1) {
            throw std::runtime_error("GuiContainer::loadFromFile: unknown element \"" + name + "\"");
        }
        //element = GuiFabric::createFromString(name);
        element = fabric[name]();
        //}
        element->loadFromFile(file);
        add(id, std::move(element));
    }
    file.read(reinterpret_cast<char*>(&m_activeObject), sizeof(m_activeObject));
    //std::shared_ptr<GUIContainer> fag;
    //fag["aefaef"];
    //fag->getContainer("agfaeg")->getContainer("aegaerg")->getInteractable("afeaef")->setSignal();
    Interactable::loadFromFile(file);
}

void GUIContainer::saveToFile (std::ofstream &file) const {
    uint64_t size;
    //size = M_NAME.size();
    //file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    //file.write(&M_NAME[0], sizeof(char) * size);
    size = m_content.size();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (const auto &el : m_content) {
        size = el.first.size();
        file.write(reinterpret_cast<const char*>(&size), sizeof(size));
        file.write(&el.first[0], sizeof(char) * size);
        el.second->saveToFile(file);
    }
    file.write(reinterpret_cast<const char*>(&m_activeObject), sizeof(m_activeObject));
    Interactable::saveToFile(file);
}