#ifndef CONSOLE_HPP
#define CONSOLE_HPP

#include <GUI/Interactable.hpp>
#include <GUI/TextBased.hpp>
#include <GUI/RenderTextureBased.hpp>
#include <GUI/InputLine.hpp>

class Console : public Interactable, public TextBased, public RenderTextureBased {
    public:
        struct Command {
            sf::String m_name, m_description;
            std::function<void (const std::vector<sf::String> &)> m_func;

            Command ();
            Command (const sf::String &name, const std::function<void (const std::vector<sf::String> &)> &func);
            ~Command ();

            void setDescription (const sf::String &description);
            void getDescription () const;

            bool operator== (const sf::String &str);
            void operator() ();
            void operator() (const std::vector<sf::String> &args);
        };
    private:
        void updateRender () override;
        void updateContent () override;
        void interact (const sf::Event &event, const sf::Vector2f &mousePos) override;
        void drawVisible (sf::RenderTarget& target, sf::RenderStates states) const override;
    public:
        Console ();
        ~Console ();

        void saveHistoryToFile (const std::string &file);
        void loadHistoryFromFile (const std::string &file);
        std::vector<sf::String> getHistory ();
        bool receivedCommand () const;
        sf::String getLastCommand ();
        void setHistoryLimit (uint64_t limit);
        void setLimitOfSymbols (uint64_t limit);
        void setUserName (const sf::String &name);
        sf::String getUserName () const;
        void addCommand (const sf::String &name, const std::function<void (const std::vector<sf::String> &)> &func);
        //void addCommand (const Command &command);
        //void setCommandList (const std::vector<Command> &list);
        void clear ();
        using GUIElement::setSize;
        //void update () override;
        void setSize (const sf::Vector2f &size) override;

        friend std::ifstream &operator>> (std::ifstream &file, Console &console);
        friend std::ofstream &operator<< (std::ofstream &file, const Console &console);
    private:
        //static const uint64_t M_HISTORY_LIMIT = 200;
        uint64_t m_lineNumber, m_limitOfSymbols, m_historyLimit;
        sf::RectangleShape m_cursor;
        std::vector<sf::String> m_history;
        std::vector<Command> m_commands;
        //std::map<sf::String, std::function<void (const std::vector<sf::String> &)>> m_commands;
        sf::String m_userName;
        sf::Clock m_clock;
        InputLine m_currentLine;
};

#endif