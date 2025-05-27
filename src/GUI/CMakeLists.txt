#add_library(${GUI_LIB} "")

include(Config.cmake)

#---------Base----------
#TextBased
add_module_lib(${TEXT_BASED_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/TextBased.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/TextBased.hpp
)
#RenderTextureBased
add_module_lib(${RENDER_TEXTURE_BASED_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/RenderTextureBased.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/RenderTextureBased.hpp
)
#GuiElement
add_module_lib(${GUI_ELEMENT_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/GUIElement.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/GUIElement.hpp
    DEPENDS
        ${SF_EXTENSIONS_LIB}
)
#Interactable
add_module_lib(${INTERACTABLE_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/Interactable.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/Interactable.hpp
    DEPENDS
        ${GUI_ELEMENT_LIB}
)
#GuiContainer
add_module_lib(${GUI_CONTAINER_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/GUIContainer.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/GUIContainer.hpp
    DEPENDS
        ${INTERACTABLE_LIB}
)
#InputLine
add_module_lib(${INPUT_LINE_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/InputLine.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/InputLine.hpp
    DEPENDS
        ${INTERACTABLE_LIB}
        ${TEXT_BASED_LIB}
)

#------Derivatives------
#Button
add_module_lib(${BUTTON_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/Button.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/Button.hpp
    DEPENDS
        ${INTERACTABLE_LIB}
        ${TEXT_BASED_LIB}
)
#Bar
add_module_lib(${BAR_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/Bar.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/Bar.hpp
    DEPENDS
        ${INTERACTABLE_LIB}
)
#Dialogue
add_module_lib(${DIALOGUE_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/Dialogue.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/Dialogue.hpp
    DEPENDS
        ${INTERACTABLE_LIB}
        ${TEXT_BASED_LIB}
)
#Slider
add_module_lib(${SLIDER_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/Slider.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/Slider.hpp
    DEPENDS
        ${INTERACTABLE_LIB}
)
#Plot2D
add_module_lib(${PLOT2D_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/Plot2D.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/Plot2D.hpp
    DEPENDS
        ${INTERACTABLE_LIB}
        ${TEXT_BASED_LIB}
        ${RENDER_TEXTURE_BASED_LIB}
        ${SF_EXTENSIONS_LIB}
        ${LEGEND_LIB}
        ${NUMERIC_METHODS_LIB}
)
#Plot3D
add_module_lib(${PLOT3D_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/Plot3D.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/Plot3D.hpp
    DEPENDS
        ${INTERACTABLE_LIB}
        ${TEXT_BASED_LIB}
        ${RENDER_TEXTURE_BASED_LIB}
)
#ScrollingWindow
add_module_lib(${SCROLLING_WINDOW_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/ScrollingWindow.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/ScrollingWindow.hpp
    DEPENDS
        ${GUI_CONTAINER_LIB}
        ${RENDER_TEXTURE_BASED_LIB}
        ${SLIDER_LIB}
)
#Selector
add_module_lib(${SELECTOR_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/Selector.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/Selector.hpp
    DEPENDS
        ${INTERACTABLE_LIB}
        ${TEXT_BASED_LIB}
        ${RENDER_TEXTURE_BASED_LIB}
)
#TextField
add_module_lib(${TEXT_FIELD_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/TextField.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/TextField.hpp
    DEPENDS
        ${INTERACTABLE_LIB}
        ${RENDER_TEXTURE_BASED_LIB}
        ${INPUT_LINE_LIB}
)
#Table
add_module_lib(${TABLE_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/Table.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/Table.hpp
    DEPENDS
        ${INTERACTABLE_LIB}
        ${TEXT_BASED_LIB}
        ${RENDER_TEXTURE_BASED_LIB}
        ${TEXT_FIELD_LIB}
)
#Page
add_module_lib(${PAGE_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/Page.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/Page.hpp
    DEPENDS
        ${GUI_CONTAINER_LIB}
)
#PageHolder
add_module_lib(${PAGE_HOLDER_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/PageHolder.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/PageHolder.hpp
    DEPENDS
        ${INTERACTABLE_LIB}
        ${PAGE_LIB}
)
#Box
add_module_lib(${BOX_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/Box.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/Box.hpp
    DEPENDS
        ${GUI_CONTAINER_LIB}
)
#Label
add_module_lib(${LABEL_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/Label.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/Label.hpp
    DEPENDS
        ${GUI_ELEMENT_LIB}
        ${TEXT_BASED_LIB}
)
#ButtonBox
add_module_lib(${BUTTON_BOX_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/ButtonBox.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/ButtonBox.hpp
    DEPENDS
        ${INTERACTABLE_LIB}
        ${RENDER_TEXTURE_BASED_LIB}
)
#Spacer
add_module_lib(${SPACER_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/Spacer.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/Spacer.hpp
    DEPENDS
        ${GUI_ELEMENT_LIB}
)
#MessageWindow
add_module_lib(${MESSAGE_WINDOW_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/MessageWindow.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/MessageWindow.hpp
    DEPENDS
        ${INTERACTABLE_LIB}
        ${TEXT_BASED_LIB}
        ${RENDER_TEXTURE_BASED_LIB}
        ${BUTTON_LIB}
)
#Cursor
add_module_lib(${CURSOR_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/Cursor.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/Cursor.hpp
    DEPENDS
        ${INTERACTABLE_LIB}
)
#RichText
add_module_lib(${RICH_TEXT_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/RichText.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/RichText.hpp
    DEPENDS
        ${GUI_ELEMENT_LIB}
)
#HeatMap
add_module_lib(${HEAT_MAP_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/HeatMap.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/HeatMap.hpp
    DEPENDS
        ${INTERACTABLE_LIB}
        ${RENDER_TEXTURE_BASED_LIB}
        ${GENERAL_LIB}
)
#Legend
add_module_lib(${LEGEND_LIB}
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/Legend.cpp
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}/Legend.hpp
    DEPENDS
        ${INTERACTABLE_LIB}
        ${RENDER_TEXTURE_BASED_LIB}
        ${TEXT_BASED_LIB}
)
#------GUI_LIB------
add_module_lib(${GUI_LIB}
    #PRIVATE
    #    ${CMAKE_CURRENT_LIST_DIR}/RichText.cpp
    #PUBLIC
    #    ${CMAKE_CURRENT_LIST_DIR}/RichText.hpp
    DEPENDS
        ${BUTTON_LIB}
        ${BAR_LIB}
        ${DIALOGUE_LIB}
        ${PLOT2D_LIB}
        ${PLOT3D_LIB}
        ${SCROLLING_WINDOW_LIB}
        ${SELECTOR_LIB}
        ${SLIDER_LIB}
        ${TEXT_FIELD_LIB}
        ${TABLE_LIB}
        ${PAGE_LIB}
        ${PAGE_HOLDER_LIB}
        ${BOX_LIB}
        ${LABEL_LIB}
        ${BUTTON_BOX_LIB}
        ${SPACER_LIB}
        ${MESSAGE_WINDOW_LIB}
        ${CURSOR_LIB}
        ${RICH_TEXT_LIB}
        ${HEAT_MAP_LIB}
)
# target_sources(${GUI_LIB}
#     PRIVATE
#         ${CMAKE_CURRENT_LIST_DIR}/TextBased.cpp
#         ${CMAKE_CURRENT_LIST_DIR}/RenderTextureBased.cpp
#         ${CMAKE_CURRENT_LIST_DIR}/Bar.cpp
#         ${CMAKE_CURRENT_LIST_DIR}/Button.cpp
#         #${CMAKE_CURRENT_LIST_DIR}/Console.cpp
#         ${CMAKE_CURRENT_LIST_DIR}/Dialogue.cpp
#         ${CMAKE_CURRENT_LIST_DIR}/Plot2D.cpp
#         ${CMAKE_CURRENT_LIST_DIR}/Plot3D.cpp
#         ${CMAKE_CURRENT_LIST_DIR}/GUIElement.cpp
#         ${CMAKE_CURRENT_LIST_DIR}/Interactable.cpp
#         ${CMAKE_CURRENT_LIST_DIR}/ScrollingWindow.cpp
#         ${CMAKE_CURRENT_LIST_DIR}/Selector.cpp
#         ${CMAKE_CURRENT_LIST_DIR}/Slider.cpp
#         ${CMAKE_CURRENT_LIST_DIR}/TextField.cpp
#         ${CMAKE_CURRENT_LIST_DIR}/Table.cpp
#         ${CMAKE_CURRENT_LIST_DIR}/InputLine.cpp
#         ${CMAKE_CURRENT_LIST_DIR}/Page.cpp
#         ${CMAKE_CURRENT_LIST_DIR}/PageHolder.cpp
#         ${CMAKE_CURRENT_LIST_DIR}/GUIContainer.cpp
#         ${CMAKE_CURRENT_LIST_DIR}/Box.cpp
#         ${CMAKE_CURRENT_LIST_DIR}/Label.cpp
#         ${CMAKE_CURRENT_LIST_DIR}/ButtonBox.cpp
#         ${CMAKE_CURRENT_LIST_DIR}/Spacer.cpp
#         ${CMAKE_CURRENT_LIST_DIR}/MessageWindow.cpp
#         ${CMAKE_CURRENT_LIST_DIR}/Cursor.cpp
#         ${CMAKE_CURRENT_LIST_DIR}/RichText.cpp
#     PUBLIC
#         ${CMAKE_CURRENT_LIST_DIR}/DrawableHandler.hpp
#         ${CMAKE_CURRENT_LIST_DIR}/TextBased.hpp
#         ${CMAKE_CURRENT_LIST_DIR}/RenderTextureBased.hpp
#         ${CMAKE_CURRENT_LIST_DIR}/Bar.hpp
#         ${CMAKE_CURRENT_LIST_DIR}/Button.hpp
#         #${CMAKE_CURRENT_LIST_DIR}/Console.hpp
#         ${CMAKE_CURRENT_LIST_DIR}/Dialogue.hpp
#         ${CMAKE_CURRENT_LIST_DIR}/Plot2D.hpp
#         ${CMAKE_CURRENT_LIST_DIR}/Plot3D.hpp
#         ${CMAKE_CURRENT_LIST_DIR}/GUIElement.hpp
#         ${CMAKE_CURRENT_LIST_DIR}/Interactable.hpp
#         ${CMAKE_CURRENT_LIST_DIR}/ScrollingWindow.hpp
#         ${CMAKE_CURRENT_LIST_DIR}/Selector.hpp
#         ${CMAKE_CURRENT_LIST_DIR}/Slider.hpp
#         ${CMAKE_CURRENT_LIST_DIR}/TextField.hpp
#         ${CMAKE_CURRENT_LIST_DIR}/Table.hpp
#         ${CMAKE_CURRENT_LIST_DIR}/InputLine.hpp
#         ${CMAKE_CURRENT_LIST_DIR}/Page.hpp
#         ${CMAKE_CURRENT_LIST_DIR}/PageHolder.hpp
#         ${CMAKE_CURRENT_LIST_DIR}/GUIContainer.hpp
#         ${CMAKE_CURRENT_LIST_DIR}/Box.hpp
#         ${CMAKE_CURRENT_LIST_DIR}/Label.hpp
#         ${CMAKE_CURRENT_LIST_DIR}/ButtonBox.hpp
#         ${CMAKE_CURRENT_LIST_DIR}/Spacer.hpp
#         ${CMAKE_CURRENT_LIST_DIR}/MessageWindow.hpp
#         ${CMAKE_CURRENT_LIST_DIR}/Cursor.hpp
#         ${CMAKE_CURRENT_LIST_DIR}/RichText.hpp
# )

# target_link_libraries(${GUI_LIB} 
#     PRIVATE
#         ${SFML_LIBS}
#         ${GENERAL_LIB}
# )

# target_include_directories(${GUI_LIB}
#     PRIVATE
#         ${MODULES_INCLUDE_PATH}
#         ${CMAKE_CURRENT_LIST_DIR}
# )

# add_definitions(${COMPILE_FLAGS})