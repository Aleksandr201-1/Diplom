#ifndef PAGE_HOLDER_HPP
#define PAGE_HOLDER_HPP

#include <GUI/Page.hpp>

class PageHolder : public Interactable {
    private:
        void updateContent () override;
        void interact (const sf::Event &event, const sf::Vector2f &mousePos) override;
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
    public:
        PageHolder ();
        PageHolder (uint64_t n);
        ~PageHolder ();

        void addPage (std::unique_ptr<Page> &&page);
        std::unique_ptr<Page> &getPage (uint64_t id);
        const std::unique_ptr<Page> &getPage (uint64_t id) const;
        void setActivePage (uint64_t id);
        void setPageCount (uint64_t count);
        void next ();
        void prev ();
        //uint64_t &getPageRef ();
        //const uint64_t &getPageRef () const;
        //void updateAll ();
        using GUIElement::setSize;
        void setSize (const sf::Vector2f &size) override;
        using GUIElement::loadFromFile;
        using GUIElement::saveToFile;
        void loadFromFile (std::ifstream &file) override;
        void saveToFile (std::ofstream &file) const override;

        std::unique_ptr<Page> &operator[] (uint64_t id);
        const std::unique_ptr<Page> &operator[] (uint64_t id) const;
    private:
        static const std::string M_NAME;
        static const GuiFabric<PageHolder> m_pageHolderFabric;
        //std::map<uint64_t, std::unique_ptr<Page>> m_pages;
        std::vector<std::unique_ptr<Page>> m_pages;
        uint64_t m_activePage;
};

#endif