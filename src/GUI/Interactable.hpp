#ifndef INTERACTABLE_HPP
#define INTERACTABLE_HPP

#include <GUI/GUIElement.hpp>
#include <General/Underlying.hpp>
#include <functional>

class Interactable : public GUIElement {
    public:
        enum Signal : uint32_t {
            PRESSED          = 1 << 0,
            RELEASED         = 1 << 1,
            TOUCHED          = 1 << 2,
            RECEIVED_FOCUS   = 1 << 3,
            LOST_FOCUS       = 1 << 4,
            RESIZED          = 1 << 5,
            RECEIVED_COMMAND = 1 << 6,
            SELECTED_TEXT    = 1 << 7,
            COPIED_TEXT      = 1 << 8,
            SELECTED_ITEM    = 1 << 9,
            CHANGED_VALUE    = 1 << 10,
            VISUAL_CHANGE    = 1 << 11,
            ERROR_RECEIVED   = 1 << 12,
            ADDED_ITEM       = 1 << 13,
            COUNT_OF_SIGNALS = 14
        };
        enum ActivationTrigger : uint8_t {
            MOUSE    = 1 << 0,
            KEYBOARB = 1 << 1,
            JOYSTICK = 1 << 2
        };
    protected:
        virtual void updateContent ();
        virtual void interact (const sf::Event &event, const sf::Vector2f &mousePos) = 0;
        void clearSignals ();
        //void addSignal (Signal signal);
        void updateFocusCondition (const sf::Event &event, const sf::Vector2f &mousePos);
        void updateSizeCondition (const sf::Event &event, const sf::Vector2f &mousePos);
    public:
        Interactable ();
        virtual ~Interactable ();
        void setInteractionFocus (bool focus);
        bool getInteractionFocus () const;
        void setSignal (Signal signal, const std::function<void (void)> &function);
        void handleSignal (Signal signal);
        void needMousePressed (bool signal);
        bool needUpdate () const;
        bool getSignalStatus (Signal signal) const;
        //const std::array<bool, static_cast<uint64_t>(Signal::COUNT_OF_SIGNALS)> &getAllSignals () const;
        uint32_t getAllSignals () const;
        void setActivationTrigger (const sf::Mouse::Button &mouse);
        void setActivationTrigger (const sf::Keyboard::Key &key);
        void disableTrigger (ActivationTrigger trigger);
        void handleEvent (const sf::Event &event, const sf::Vector2f &mousePos);
        void update ();
        void loadFromFile (std::ifstream &file) override;
        void saveToFile (std::ofstream &file) const override;
    protected:
        uint8_t m_trigger;
        sf::Mouse::Button m_mouse;
        sf::Keyboard::Key m_key;
        bool m_interactFocus, m_resizeFocus;
        bool m_pressInsideToLoseFocus, m_pressOutsideToLoseFocus;
        bool m_needTriggerPressed, m_needTriggerInside;
        bool m_needUpdate;
        uint32_t m_signals;
        //std::array<bool, static_cast<uint64_t>(Signal::COUNT_OF_SIGNALS)> m_signals; //need change to 4 or 8 byte value lol
        std::map<Signal, std::function<void (void)>> m_workers;
};

#endif