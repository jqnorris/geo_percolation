TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11

SOURCES += main.cpp \
    Bond.cpp \
    Site.cpp \
    Strength.cpp \
    Lattice.cpp \
    Algorithm.cpp

include(deployment.pri)
qtcAddDeployment()

HEADERS += \
    Abstract_Classes.h \
    Forward_Declarations.h

