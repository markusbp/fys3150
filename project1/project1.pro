TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        ../../ComputationalPhysics/doc/Programs/LecturePrograms/programs/cppLibrary/lib.cpp \
        main.cpp

LIBS += -larmadillo

HEADERS +=
        ../../ComputationalPhysics/doc/Programs/LecturePrograms/programs/cppLibrary/lib.h \
