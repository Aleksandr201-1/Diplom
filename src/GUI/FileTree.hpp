#ifndef FILE_TREE_HPP
#define FILE_TREE_HPP

#include <GUI/Interactable.hpp>
#include <GUI/TextBased.hpp>
#include <GUI/RenderTextureBased.hpp>
#include <memory>
#include <list>

class FileTree : public Interactable, public TextBased, public RenderTextureBased {
    private:
        class FileNode {
            //private:
                //void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
            public:
                FileNode ();
                ~FileNode ();
            private:
                //static float m_Hoffset, m_Woffset;
                //static std::map<sf::String, FileNode> m_childs;
                sf::String m_name;
                uint64_t m_level;
                bool m_drawChilds;
                std::vector<FileNode> m_childs;
                FileNode *m_parent;
                //std::vector<std::shared_ptr<FileNode>> m_childs;
        };
    private:
        void updateRender () override;
        void updateContent () override;
        void interact (const sf::Event &event, const sf::Vector2f &mousePos) override;
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
    public:
        FileTree ();
        ~FileTree ();

        void addItem (const sf::String &str, uint64_t level, uint64_t pos);
        void addAfterActive (const sf::String &str);
        void deleteItem ();
        void addNewSubGroup ();
        uint64_t getCurrentLevel () const;
        uint64_t getCurrentId () const;

        using GUIElement::setSize;
        void setSize (const sf::Vector2f &size) override;
    private:
        float m_Hoffset, m_Woffset;
        uint64_t m_currId;
        std::vector<FileNode> m_nodes;
        FileNode *m_parent;
        //std::list<FileNode> m_nodes;
        //std::shared_ptr<FileNode> m_root;
};

#endif