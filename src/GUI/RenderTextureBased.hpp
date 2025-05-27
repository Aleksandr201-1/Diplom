#ifndef RENDER_TEXTURE_BASED_HPP
#define RENDER_TEXTURE_BASED_HPP

#include <SFML/Graphics.hpp>
#include <fstream>

class RenderTextureBased {
    protected:
        virtual void recalculateSize () = 0;
        virtual void updateRender () = 0;
    public:
        RenderTextureBased ();
        RenderTextureBased (uint64_t width, uint64_t height);
        RenderTextureBased (const sf::Vector2u &size);
        virtual ~RenderTextureBased ();

        void setRenderSize (uint64_t width, uint64_t height);
        void setRenderSize (const sf::Vector2u &size);
        void setFonColor (const sf::Color &color);

        void loadFromFile (std::ifstream &file);
        void saveToFile (std::ofstream &file) const;
    protected:
        sf::Color m_clearColor;
        bool m_needRedraw;
        sf::RenderTexture m_render;
        sf::Sprite m_visibleContent;
};

#endif