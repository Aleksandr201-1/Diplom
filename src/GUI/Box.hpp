#ifndef BOX
#define BOX

#include <GUI/GUIContainer.hpp>
#include <memory>
#include <map>

class Box : public GUIContainer {
    public:
        enum class Orientation {
            VERTICAL,
            HORIZONTAL
        };
        enum class Alignment {
            UP,
            DOWN,
            LEFT,
            RIGHT,
            CENTER
        };
    private:
        void recalculateSize () override;
        void realign ();
        void updateContent () override;
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
    public:
        Box ();
        Box (uint64_t offset, Alignment alignment, Orientation orientation);
        ~Box ();

        void setOffset (uint64_t offset);
        void setAlignment (Alignment alignment);
        void setOrientation (Orientation orientation);
        void needUniteSize (bool status);
        void saveOffset (bool status);
        //void addGUIElement (std::shared_ptr<GUIElement> &&el);
        //bool needUpdate () const;
        //void update () override;
        using GUIElement::setSize;
        void setSize (const sf::Vector2f &size) override;
        void loadFromFile (std::ifstream &file) override;
        void saveToFile (std::ofstream &file) const override;
    private:
        static const std::string M_NAME;
        static const GuiFabric<Box> m_boxFabric;
        const uint64_t M_MINIMAL_OFFSET = 5;
        float m_maxWidth, m_maxHeight;
        uint64_t m_offset;
        Alignment m_alignment;
        Orientation m_orientation;
        bool m_needRealign, m_needUniteSize, m_saveOffset;
        //sf::RectangleShape m_border;
};

#endif