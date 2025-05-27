#ifndef TABLE_HPP
#define TABLE_HPP

#include <GUI/Interactable.hpp>
#include <GUI/TextBased.hpp>
#include <GUI/RenderTextureBased.hpp>
#include <GUI/TextField.hpp>
#include <cmath>

class Table : public Interactable, public TextBased, public RenderTextureBased {
    private:
        void recalculateSize () override;
        void updateRender () override;
        void updateContent () override;
        void interact (const sf::Event &event, const sf::Vector2f &mousePos) override;
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
    public:
        Table ();
        Table (const std::vector<std::vector<sf::String>> &content);
        ~Table ();

        void setContent (const std::vector<std::vector<sf::String>> &content);
        void setCols (uint64_t cols);
        void setRows (uint64_t rows);
        void addRow (uint64_t i, const std::vector<sf::String> &row);
        uint64_t getRowCount () const;
        uint64_t getColCount () const;
        void setColSize (uint64_t i, float size);
        void setRowSize (uint64_t i, float size);
        void changeItem (uint64_t i, uint64_t j, const sf::String &item);
        void setGrid (bool grid);
        void editStatus (bool status);
        //void update () override;
        using GUIElement::setSize;
        void setSize (const sf::Vector2f &size) override;
        void loadFromFile (std::ifstream &file) override;
        void saveToFile (std::ofstream &file) const override;

        //std::vector<sf::String> &operator[] (uint64_t idx);
        const std::vector<sf::String> &operator[] (uint64_t i) const;
    private:
        static const std::string M_NAME;
        static const GuiFabric<Table> m_tableFabric;
        std::vector<std::vector<sf::String>> m_content;
        std::vector<float> m_HLengths, m_VLengths;
        bool m_grid, m_editable;
        TextField m_field;
        sf::RectangleShape m_shape;
};

#endif