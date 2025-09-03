# Nombre del módulo
TARGET = integration

# Compilador y flags
CXX = c++
CXXFLAGS = -O3 -Wall -shared -std=c++17 -fPIC

# Includes de Python y pybind11
PYTHON_INCLUDES = $(shell python3 -m pybind11 --includes)
PYTHON_EXT = $(shell python3-config --extension-suffix)

# Archivos fuente
SRC = integration.cpp

# Regla por defecto
all: $(TARGET)$(PYTHON_EXT)

$(TARGET)$(PYTHON_EXT): $(SRC)
	$(CXX) $(CXXFLAGS) $(PYTHON_INCLUDES) $(SRC) -o $@

clean:
	rm -f $(TARGET)*.so $(TARGET)*.pyd
