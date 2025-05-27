#ifndef LEGEND_HPP
#define LEGEND_HPP

#include <GUI/Interactable.hpp>
#include <GUI/TextBased.hpp>
#include <GUI/RenderTextureBased.hpp>
#include <SFExtensions/Line.hpp>

class Legend : public Interactable, public TextBased, public RenderTextureBased {
    public:
        struct Record {
            Line line;
            sf::String string;

            Record ();
            Record (const Line &line, const sf::String &string);
            ~Record ();
        };
    private:
        void recalculateSize () override;
        void updateRender () override;
        void updateContent () override;
        void interact (const sf::Event &event, const sf::Vector2f &mousePos) override;
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
    public:
        Legend ();
        ~Legend ();
        void add (const Record &record);
        void clear ();

        Record &operator[] (uint64_t i);
        Legend &operator= (const Legend &legend) = default;
        Legend &operator= (Legend &&legend) = default;
        using GUIElement::setSize;
        void setSize (const sf::Vector2f &size) override;
        void loadFromFile (std::ifstream &file) override;
        void saveToFile (std::ofstream &file) const override;
    private:
        static const std::string M_NAME;
        static const GuiFabric<Legend> m_legendFabric;
        std::vector<Record> m_records;
};

class LegendAvailable {
    public:
        LegendAvailable ();
        virtual ~LegendAvailable ();
        virtual std::shared_ptr<Legend> generateLegend () = 0;
};

#endif