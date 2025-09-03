#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <cmath>
#include <iostream>

namespace py = pybind11;
// ------------------
// 1) Trapecio con función Python (menos eficiente)
// ------------------
double trapezoid(py::function fx, double x0, double xN, int N) {
    double h = (xN - x0) / N;
    double sum = 0.0;

    sum += 0.5 * (fx(x0).cast<double>() + fx(xN).cast<double>());
    for (int i = 1; i < N; i++) {
        double xi = x0 + i * h;
        sum += fx(xi).cast<double>();
    }
    return sum * h;
}

// ------------------
// 2) Trapecio con arrays de NumPy (más eficiente)
// ------------------
double trapezoid_array(py::array_t<double> x, py::array_t<double> y) {
    auto buf_x = x.unchecked<1>();
    auto buf_y = y.unchecked<1>();
    int N = buf_x.shape(0);

    if (N != buf_y.shape(0)) {
        throw std::runtime_error("x and y must have the same length");
    }

    double h = (buf_x(N-1) - buf_x(0)) / (N-1);
    double sum = 0.5 * (buf_y(0) + buf_y(N-1));

    for (int i = 1; i < N-1; i++) {
        sum += buf_y(i);
    }

    return sum * h;
}

double simpsons_array(py::array_t<double> x, py::array_t<double> y) {
    auto buf_x = x.unchecked<1>();
    auto buf_y = y.unchecked<1>();
    int N = buf_x.shape(0);

    if (N != buf_y.shape(0)) {
        throw std::runtime_error("x and y must have the same length");
    }
    if (N % 2 == 0) {
        throw std::runtime_error("N must be odd for Simpson's rule (even number of intervals).");
    }

    double h = (buf_x[N-1] - buf_x[0]) / (N - 1);

    double sum = buf_y[0] + buf_y[N-1];
    for (int i = 1; i < N-1; i++) {
        if (i % 2 == 0) {  // even index
            sum += 2 * buf_y[i];
        } else {           // odd index
            sum += 4 * buf_y[i];
        }
    }

    return sum * h / 3.0;
}


PYBIND11_MODULE(integration, m) {
    m.doc() = "Módulo de integración numérica con regla del trapecio";
    m.def("trapezoid", &trapezoid,
          R"pbdoc(
          Aproxima la integral definida de una función en [x0, xN]
          usando la regla del trapecio.
          Parámetros
          ----------
          fx : callable
              Función f(x).
          x0 : float
              Límite inferior.
          xN : float
              Límite superior.
          N : int
              Número de subdivisiones.
          Retorna
          -------
          float
              Aproximación de la integral.
          )pbdoc",
          py::arg("fx"), py::arg("x0"), py::arg("xN"), py::arg("N"));
    m.def("trapezoid_array", &trapezoid_array,
          R"pbdoc(
          Aproxima la integral definida a partir de valores ya evaluados.
          Parámetros
          ----------
          x : numpy.ndarray
              Puntos en el eje x (1D, ordenados).
          y : numpy.ndarray
              Valores de f(x) en esos puntos.
          Retorna
          -------
          float
              Aproximación de la integral.
          )pbdoc",
          py::arg("x"), py::arg("y"));
    m.def("Simpson_array",&simpsons_array,
    R"pbdoc(
          Aproxima la integral definida a partir de valores ya evaluados.
          Parámetros
          ----------
          x : numpy.ndarray
              Puntos en el eje x (1D, ordenados).
          y : numpy.ndarray
              Valores de f(x) en esos puntos.
          Retorna
          -------
          float
              Aproximación de la integral.
          )pbdoc",
        py::arg("x"), py::arg("y"));
}
