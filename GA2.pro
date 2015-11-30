HEADERS       = glwidget.h \
    plot.h \
    GeneticAlgorithm.h
SOURCES       = glwidget.cpp \
                main.cpp \
    plot.cpp \
    GeneticAlgorithm.cpp

QT           += widgets

# install
target.path = $$[QT_INSTALL_EXAMPLES]/opengl/hellogl2
INSTALLS += target
