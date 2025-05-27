#include <SFExtensions/SpriteManager.hpp>

SpriteManager::SpriteManager (const std::string &path, uint64_t col, uint64_t row, const sf::Color &mask) : col(col), row(row) {
    loadTexture(path, mask);
}

SpriteManager::~SpriteManager () {}

void SpriteManager::loadTexture (const std::string &path, const sf::Color &mask) {
    image.loadFromFile(path);
    image.createMaskFromColor(mask);
    texture.loadFromImage(image);
    uint64_t textureW = texture.getSize().x, textureH = texture.getSize().y;
    spriteW = textureW / col, spriteH = textureH / row;
    sprites.clear();
    for (uint64_t i = 0; i < col * row; ++i) {
        sprites.push_back(sf::Sprite(texture, sf::IntRect((i % col) * spriteW, (i / col) * spriteH, spriteW, spriteH)));
    }
}

uint64_t SpriteManager::cols () const {
    return col;
}

uint64_t SpriteManager::rows () const {
    return row;
}

sf::Sprite SpriteManager::getCustomSprite (uint64_t x, uint64_t y, uint64_t w, uint64_t h) const {
    return sf::Sprite(texture, sf::IntRect(x, y, w, h));
}

// sf::Sprite SpriteManager::operator[] (uint64_t i) const {
//     if (i > col * row) {
//         throw std::logic_error("Texture array: out of range");
//     }
//     return sf::Sprite(texture, sf::IntRect((i % col) * spriteW, (i / col) * spriteH, spriteW, spriteH));
// }

const sf::Sprite &SpriteManager::operator[] (uint64_t i) const {
    if (i > col * row) {
        throw std::logic_error("Texture array: out of range");
    }
    return sprites[i];
}

sf::Sprite SpriteManager::at (uint64_t i, uint64_t j) const {
    if (i > col || j > row) {
        throw std::logic_error("Texture array: out of range");
    }
    return sf::Sprite(texture, sf::IntRect(i * spriteW, j * spriteH, spriteW, spriteH));
}

sf::Sprite SpriteManager::at (uint64_t i) const {
    if (i > col * row) {
        throw std::logic_error("Texture array: out of range");
    }
    return sf::Sprite(texture, sf::IntRect((i % col) * spriteW, (i / col) * spriteH, spriteW, spriteH));
}