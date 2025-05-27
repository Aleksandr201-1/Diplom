//  key_id: YCAJEwfrPLaSTDKiOZr_i4a4V
//  secret: YCOZfsA-EtG0x086rANZ2hIwwruR-NeL1lcKRVc1

#ifndef LABEL_HPP
#define LABEL_HPP

#include <GUI/GUIElement.hpp>
#include <GUI/TextBased.hpp>

class Label : public GUIElement, public TextBased {
    private:
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
    public:
        Label ();
        Label (const sf::String &string);
        ~Label ();
        //sf::Uint32 &operator[] (uint64_t i);
        sf::Uint32 operator[] (uint64_t i) const;

        using GUIElement::setSize;
        void setSize (const sf::Vector2f &size) override;
        void loadFromFile (std::ifstream &file) override;
        void saveToFile (std::ofstream &file) const override;
    private:
        static const std::string M_NAME;
        static const GuiFabric<Label> m_labelFabric;
};

#endif