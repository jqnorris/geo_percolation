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
    Algorithm.cpp \
    Analysis.cpp \
    Tools.cpp \
    Simulation.cpp \
    Common_Runs.cpp

include(deployment.pri)
qtcAddDeployment()

HEADERS += \
    Abstract_Classes.h \
    Forward_Declarations.h \
    Tools.h \
    Simulation.h \
    Common_Runs.h \
    Analysis.h

