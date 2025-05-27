#include <GUI/RichText.hpp>

const std::string RichText::M_NAME = "RichText";
const RichText::GuiFabric<RichText> RichText::m_richTextFabric(RichText::M_NAME);

RichText::Style::Style ()
    : font(nullptr), color(sf::Color::Black), fontSize(28), fontStyle(sf::Text::Style::Regular) {}

RichText::Style::Style (const sf::Font &font, const sf::Color color, uint64_t fontSize, uint32_t fontStyle)
    : font(&font), color(color), fontSize(fontSize), fontStyle(fontStyle) {}

RichText::Style::Style (const Style &style)
    : font(style.font), color(style.color), fontSize(style.fontSize), fontStyle(style.fontStyle) {}

RichText::Style::Style (Style &&style)
    : font(style.font), color(style.color), fontSize(style.fontSize), fontStyle(style.fontStyle) {}

RichText::Style::~Style () {}

RichText::Style &RichText::Style::operator= (const Style &s) {
    font = s.font;
    color = s.color;
    fontSize = s.fontSize;
    fontStyle = s.fontStyle;
    return *this;
}

RichText::Style &RichText::Style::operator= (Style &&s) {
    font = s.font;
    color = s.color;
    fontSize = s.fontSize;
    fontStyle = s.fontStyle;
    return *this;
}

bool operator== (const RichText::Style &s1, const RichText::Style &s2) {
    return s1.color == s2.color &&
           s1.fontSize == s2.fontSize &&
           s1.fontStyle == s2.fontStyle &&
           s1.font->getInfo().family == s2.font->getInfo().family;
}

bool operator!= (const RichText::Style &s1, const RichText::Style &s2) {
    return !(s1 == s2);
}

RichText::StyleInterval::StyleInterval () {}

RichText::StyleInterval::StyleInterval (const Style &style, uint64_t begin, uint64_t end)
    : style(style), begin(begin), end(end) {}

RichText::StyleInterval::StyleInterval (const StyleInterval &styleInterval)
    : style(styleInterval.style), begin(styleInterval.begin), end(styleInterval.end) {}

RichText::StyleInterval::StyleInterval (StyleInterval &&styleInterval)
    : style(std::move(styleInterval.style)), begin(styleInterval.begin), end(styleInterval.end) {}

RichText::StyleInterval::~StyleInterval () {}

RichText::StyleInterval &RichText::StyleInterval::operator= (const RichText::StyleInterval &s) {
    style = s.style;
    begin = s.begin;
    end = s.end;
    return *this;
}

RichText::StyleInterval &RichText::StyleInterval::operator= (StyleInterval &&s) {
    style = s.style;
    begin = s.begin;
    end = s.end;
    return *this;
}

bool operator== (const RichText::StyleInterval &s1, const RichText::StyleInterval &s2) {
    return s1.style == s2.style &&
           s1.begin == s2.begin &&
           s1.end == s2.end;
}

bool operator!= (const RichText::StyleInterval &s1, const RichText::StyleInterval &s2) {
    return !(s1 == s2);
}

void RichText::calcMaxHeight () {
    uint64_t height = 0;
    for (uint64_t i = 0; i < m_intervals.size(); ++i) {
        height = std::max(height, m_intervals[i].style.fontSize);
    }
    m_maxHeight = height;
}

void RichText::drawVisible (sf::RenderTarget& target, sf::RenderStates states) const {
    sf::Vertex lines[2];
    lines[0].color = lines[1].color = sf::Color::Black;
    lines[0].position.y = lines[1].position.y = m_maxHeight + m_lineOffset;
    //lines[0].position.x = 0;
    float width = 0;
    float verticalOffset = 0;
    bool underlined = false;
    //for (const auto &interval : m_intervals) {
    for (uint64_t i = 0; i < m_intervals.size(); ++i) {
        const auto &interval = m_intervals[i];
        verticalOffset = m_maxHeight - interval.style.fontSize;
        
        m_text.setFont(*interval.style.font);
        m_text.setFillColor(interval.style.color);
        m_text.setCharacterSize(interval.style.fontSize);
        if (interval.style.fontStyle & sf::Text::Style::Underlined) {
            underlined = true;
            m_text.setStyle(interval.style.fontStyle - sf::Text::Style::Underlined);
        } else {
            m_text.setStyle(interval.style.fontStyle);
        }
        m_text.setString(m_string.substring(interval.begin, interval.end - interval.begin + 1));
        m_text.setPosition(width, verticalOffset);
        target.draw(m_text, states);
        lines[0].position.x = width;
        width += m_text.getGlobalBounds().width;
        if (i != m_intervals.size() - 1) {
            uint64_t idx1, idx2;
            idx1 = m_intervals[i].end;
            idx2 = m_intervals[i + 1].begin;
            float candidate1 = m_intervals[i].style.font->getKerning(m_string[idx1], m_string[idx2], m_intervals[i].style.fontSize);
            float candidate2 = m_intervals[i + 1].style.font->getKerning(m_string[idx1], m_string[idx2], m_intervals[i + 1].style.fontSize);
            width += std::max(candidate1, candidate2);
        }
        lines[1].position.x = width;
        lines[0].color = lines[1].color = interval.style.color;
        if (underlined) {
            for (uint64_t i = 0; i < m_lineThickness; ++i) {
                target.draw(lines, 2, sf::Lines, states);
                lines[0].position.x += 1.f;
                lines[1].position.x += 1.f;
            }
            underlined = false;
        }
    }
}

RichText::RichText () {
    m_intervals.push_back(StyleInterval());
    m_needBorders = false;
}

RichText::RichText (const sf::String &string, const RichText::Style &style) {
    m_intervals.push_back(StyleInterval(style, 0, string.getSize() - 1));
    m_string = string;
    m_size = {100, 28};
    m_maxHeight = 28.f;
    m_lineOffset = 2.f;
    m_lineThickness = 2;
}

RichText::~RichText () {

}

void RichText::setFont (const sf::Font &font, uint64_t begin, uint64_t end) {

}

void RichText::setColor (const sf::Color &color, uint64_t begin, uint64_t end) {

}

void RichText::setFontSize (uint64_t fontSize, uint64_t begin, uint64_t end) {

}

void RichText::setFontStyle (uint32_t fontStyle, uint64_t begin, uint64_t end) {

}

void RichText::setStyleInterval (Style style, uint64_t begin, uint64_t end) {
    if (begin == end) {
        return;
    }
    if (begin > end) {
        throw std::logic_error("RichText::setStyleInterval: begin pos cant be greater than end");
    }
    std::vector<StyleInterval> newIntervals;
    for (uint64_t i = 0; i < m_intervals.size(); ++i) {
        StyleInterval &el = m_intervals[i];
        if (el.end <= begin) {
            newIntervals.push_back(el);
            continue;
        }
        if (el.end > begin && el.end <= end) {
            if (el.begin < begin) {
                el.end = begin;
                ++begin;
                newIntervals.push_back(el);
            }
            continue;
        }
        if (el.begin == begin && el.end < end) {
            continue;
        }
        if (el.end > end) {
            
        }
        if (el.begin <= begin && el.end >= end) {
            StyleInterval &prev = el, next(el);
            prev.end = begin - 1;
            next.begin = end + 1;
            m_intervals.insert(m_intervals.begin() + i + 1, {StyleInterval(style, begin, end), next});
            i += 2;
        } else if (el.begin > begin && el.end < end) {

        } else if (el.end > begin && el.begin < begin) {

        }
    }
    m_intervals = newIntervals;
    if (style.fontSize > m_maxHeight) {
        m_maxHeight = style.fontSize;
    }
    calcMaxHeight();
}

void RichText::setUnderlineThickness (uint64_t thickness) {
    m_lineThickness = thickness;
}

void RichText::setUnderlineOffset (float offset) {
    m_lineOffset = offset;
}

void RichText::setSize (const sf::Vector2f &size) {
    m_size = size;
}

void RichText::loadFromFile (std::ifstream &file) {
    GUIElement::loadFromFile(file);
}

void RichText::saveToFile (std::ofstream &file) const {
    GUIElement::saveToFile(file);
}