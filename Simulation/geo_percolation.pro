TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

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

