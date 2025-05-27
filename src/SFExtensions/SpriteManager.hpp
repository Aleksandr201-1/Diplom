#ifndef SPRITE_MANAGER
#define SPRITE_MANAGER

#include <SFML/Graphics.hpp>

class SpriteManager {
    public:
        SpriteManager (const std::string &path, uint64_t col, uint64_t row, const sf::Color &mask = sf::Color::Transparent); //cols - столбцы, rows - строки
        ~SpriteManager ();
        void loadTexture (const std::string &path, const sf::Color &mask = sf::Color::Transparent);
        uint64_t cols () const;
        uint64_t rows () const;
        sf::Sprite getCustomSprite (uint64_t x, uint64_t y, uint64_t w, uint64_t h) const;
        //sf::Sprite operator[] (uint64_t i);
        const sf::Sprite &operator[] (uint64_t i) const;
        sf::Sprite at (uint64_t i, uint64_t j) const;
        sf::Sprite at (uint64_t i) const;
    private:
        sf::Image image;
        sf::Texture texture;
        std::vector<sf::Sprite> sprites;
        uint64_t spriteW, spriteH;
        uint64_t col, row;
};

#endif