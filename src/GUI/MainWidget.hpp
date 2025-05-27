#ifndef MAIN_WIDGET_HPP
#define MAIN_WIDGET_HPP

#include <GUI/PageHolder.hpp>
#include <GUI/Cursor.hpp>
#include <GUI/MessageWindow.hpp>

class MainWidget {
    private:
    public:
        MainWidget ();
        ~MainWidget ();

        void run ();
    private:
        PageHolder m_pageHolder;
        Cursor m_cursor;
        MessageWindow m_message;
        sf::RenderWindow window;
};

#endif