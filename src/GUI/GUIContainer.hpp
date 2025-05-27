#ifndef GUI_CONTAINER
#define GUI_CONTAINER

#include <GUI/Interactable.hpp>
//#include <GUI/GuiFabric.hpp>
#include <stdexcept>
#include <memory>
#include <unordered_map>

//template <class T>
class GUIContainer : public Interactable {
    protected:
        virtual void recalculateSize () = 0;
        void updateContent () override;
        void interact (const sf::Event &event, const sf::Vector2f &mousePos) override;
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
    public:
        GUIContainer ();
        virtual ~GUIContainer ();

        void add (const std::string &name, const std::string &element);
        void add (const std::string &name, std::shared_ptr<GUIElement> &&element);
        template <class T, class... Args>
        std::shared_ptr<T> add (const std::string &name, Args &&...args) {
            std::shared_ptr<T> el = std::make_shared<T>(std::forward<Args>(args)...), ret = el;
            add(name, std::move(el));
            return ret;
        }
        template <class T>
        std::shared_ptr<T> get (const std::string &name) {
            return std::dynamic_pointer_cast<T>(m_content.at(name));
        }
        void remove (const std::string &name);
        std::shared_ptr<GUIElement> &getElement (const std::string &name);
        const std::shared_ptr<GUIElement> &getElement (const std::string &name) const;
        std::shared_ptr<GUIContainer> getContainer (const std::string &name);
        std::shared_ptr<const GUIContainer> getContainer (const std::string &name) const;
        std::shared_ptr<Interactable> getInteractable (const std::string &name);
        std::shared_ptr<const Interactable> getInteractable (const std::string &name) const;
        uint64_t getContentSize () const;
        void forceRecalculateSize ();
        std::shared_ptr<GUIContainer> at (const std::string &name);
        std::shared_ptr<const GUIContainer> at (const std::string &name) const;
        //void update () override;
        //void setSize (const sf::Vector2f &size);

        std::shared_ptr<GUIContainer> operator[] (const std::string &name);
        std::shared_ptr<const GUIContainer> operator[] (const std::string &name) const;
        // std::shared_ptr<GUIContainer> operator[] (uint64_t i);
        // const std::shared_ptr<GUIContainer> operator[] (uint64_t i) const;
        //std::shared_ptr<GUIContainer> &operator() (const std::string &name);
        //const std::shared_ptr<GUIContainer> &operator() (const std::string &name) const;
        void loadFromFile (std::ifstream &file) override;
        void saveToFile (std::ofstream &file) const override;
    protected:
        // static const std::string M_NAME;
        // static const uint64_t M_ITEM_LIMIT = 10;
        // std::array<sf::String, std::shared_ptr<GUIElement>, M_ITEM_LIMIT> m_elements;
        // static std::map<sf::String, std::shared_ptr<GUIElement>> m_content;
        std::unordered_map<std::string, std::shared_ptr<GUIElement>> m_content;
        std::vector<std::string> m_names;
        uint64_t m_activeObject;
        //bool m_needRecalculateSize;
        //std::vector<std::pair<std::string, std::shared_ptr<GUIElement>>> m_content;
        //uint64_t m_activeObject;
};

#endif