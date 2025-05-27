#include <PageHolder.hpp>

const std::string PageHolder::M_NAME = "PageHolder";
const PageHolder::GuiFabric<PageHolder> PageHolder::m_pageHolderFabric(PageHolder::M_NAME);

void PageHolder::updateContent () {
    if (!m_pages.empty()) {
        m_pages[m_activePage]->update();
    }
}

void PageHolder::interact (const sf::Event &event, const sf::Vector2f &mousePos) {
    if (!m_pages.empty()) {
        m_pages[m_activePage]->handleEvent(event, mousePos);
    }
    //m_needUpdate = m_needUpdate || m_pages[m_activePage]->needUpdate();
}

void PageHolder::drawVisible (sf::RenderTarget& target, sf::RenderStates states) const {
    if (!m_pages.empty()) {
        target.draw(*m_pages.at(m_activePage), states);
    }
}

PageHolder::PageHolder () {
    m_activePage = 0;
}

PageHolder::PageHolder (uint64_t n) : PageHolder() {
    m_pages.resize(n);
    for (uint64_t i = 0; i < n; ++i) {
        m_pages[i] = std::make_unique<Page>();
    }
}

PageHolder::~PageHolder () {}

void PageHolder::addPage (std::unique_ptr<Page> &&page) {
    //m_pages.insert(std::make_pair(id, std::move(page)));
    m_pages.push_back(std::move(page));
    //m_pages[m_pages.size()] = std::move(page);
}

std::unique_ptr<Page> &PageHolder::getPage (uint64_t id) {
    return m_pages[id];
}

const std::unique_ptr<Page> &PageHolder::getPage (uint64_t id) const {
    return m_pages.at(id);
}

void PageHolder::setActivePage (uint64_t id) {
    m_activePage = id;
}

void PageHolder::setPageCount (uint64_t count) {
    uint64_t oldCount = m_pages.size();
    m_pages.resize(count);
    for (uint64_t i = oldCount; i < count; ++i) {
        m_pages[i] = std::make_unique<Page>();
    }
}

void PageHolder::next () {
    m_activePage = m_pages.size() - 1 > m_activePage ? m_activePage + 1 : 0;
}

void PageHolder::prev () {
    m_activePage = m_activePage > 0 ? m_activePage - 1 : m_pages.size() - 1;
}

// uint64_t &PageHolder::getPageRef () {
//     return m_activePage;
// }

// const uint64_t &PageHolder::getPageRef () const {
//     return m_activePage;
// }

// void PageHolder::updateAll () {
//     for (auto &el : m_pages) {
//         el->update();
//         //el.second->update();
//     }
// }

void PageHolder::setSize (const sf::Vector2f &size) {
    m_size = size;
}

void PageHolder::loadFromFile (std::ifstream &file) {
    uint64_t size;
    file.read(reinterpret_cast<char*>(&size), sizeof(size));
    //m_pages.resize(size);
    for (uint64_t i = 0; i < size; ++i) {
        auto page = std::make_unique<Page>();
        page->loadFromFile(file);
        m_pages.push_back(std::move(page));
        //m_pages[i] = std::make_unique<Page>();
        //m_pages[i]->loadFromFile(file);
    }
    Interactable::loadFromFile(file);
}

void PageHolder::saveToFile (std::ofstream &file) const {
    uint64_t size = m_pages.size();;
    //size = M_NAME.size();
    //file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    //file.write(&M_NAME[0], sizeof(char) * size);
    size = m_pages.size();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));
    for (uint64_t i = 0; i < size; ++i) {
        m_pages[i]->saveToFile(file);
    }
    Interactable::saveToFile(file);
}

std::unique_ptr<Page> &PageHolder::operator[] (uint64_t id) {
    return m_pages[id];
}

const std::unique_ptr<Page> &PageHolder::operator[] (uint64_t id) const {
    return m_pages.at(id);
}