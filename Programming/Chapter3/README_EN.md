# Spline Project Overview

Recommand to use multi-threaded compilation with `make -j4`, or it may be slow.
This project implements various spline interpolation and curve fitting methods. The project can be built, executed, and documented easily using the `make -j4` command. The components of `make` are:

```bash
make all && make run && make plot && make doxygen
```

## File Structure

### Source Files
- **`src/`**: Contains all the source code and header files.
  - **`include/`**: Directory for all header files.
    - `PPoly.hpp`: Class for piecewise polynomials.
    - `PPoly.tpp`: Template implementation for the piecewise polynomial class.
    - `PPInterpolate.hpp`: Class for PP-form spline interpolation.
    - `PPInterpolate.tpp`: Template implementation for PP-form spline interpolation.
    - `BSpline.hpp`: Class for B-splines.
    - `BSpline.tpp`: Template implementation for B-splines.
    - `BSInterpolate.hpp`: Class for B-spline interpolation.
    - `BSInterpolate.tpp`: Template implementation for B-spline interpolation.
    - `Curve.hpp`: Class for curve fitting.
    - `Curve.tpp`: Template implementation for curve fitting.
    - `BallFunction.hpp`: Class for spherical functions, including stereographic projection and spherical coordinates.
    - `MinLeastSquare.hpp`: Linear least squares class, used for convergence analysis.
  - **`src/`**: Implementation files for template classes and specific problems.
    - `PPInterpolate.cc`: Instantiations for 1st, 2nd, and 3rd-order PP spline interpolation.
    - `BSInterpolate.cc`: Instantiations for 1st, 2nd, and 3rd-order B-spline interpolation.
    - `Curve.cc`: Instantiations for 2nd and 3rd-order curve fitting.
    - `A.cc`: Implementation of Problem A.
    - `C.cc`: Implementation of Problem C.
    - `D.cc`: Implementation of Problem D.
    - `E.cc`: Implementation of Problem E.
    - `F.cc`: Implementation of Problem F.
    - `test.cc`: Tests for high-order splines.
    - `HighOrder.cc`: Additional high-order spline tests.
    - `OrderAnalysis.cc`: Code for analyzing convergence order.
    - `SplineTest.cc`: Additional spline test cases.
    - `testjson.cc`: Code to test `jsoncpp` functionality.
    - `control.json`: Input file for `testjson.cc`.
    - `plot.py`: Python script for plotting results.

### Other Files
- **`figure/`**: Directory for storing generated figures.
- **`doc/`**:
  - `design_html/`: HTML documentation generated by `doxygen`.
    - `index.html`: Main page, includes all mathematical formulas.
  - `design_pdf/`: PDF documentation generated by `doxygen`.
    - `design.pdf`: Full documentation with mathematical formulas.
  - `report.tex`: Project report in LaTeX.
  - `report.pdf`: Compiled project report in PDF format.
  - `Makefile`: Makefile for generating the report.
- **`bin/`**: Directory for storing all binary files.
- **`Makefile`**: Makefile for building and running the project.
- **`README.md`**: Project README file.
- **`Doxyfile`**: Configuration file for `doxygen`.
- **`mainpage.dox`**: Main page for `doxygen` documentation.

## Usage

The general usage framework is illustrated below and further explained in the `doxygen`-generated documentation.

```cpp
std::vector<double> x = {0, 1, 2, 3, 4, 5};
std::vector<double> y = {0, 1, 2, 3, 4, 5};
std::vector<double> boundary_condition = {0, 0};

PPInterpolate<3, double> pp(x, y, 0, boundary_condition);
BInterpolate<3, double> bs(x, y, 0, boundary_condition);
```

### Explanation of Parameters

- **`PPInterpolate`** and **`BInterpolate`**:
  - The first template parameter is the spline order (e.g., 1st, 2nd, or 3rd).
  - The second template parameter is the data type (default: `double`). You can also use `long double`, `mpf_class`, etc.

- **Constructor Parameters**:
  - **`method`**: Defines the boundary conditions.
    - For $N=1$, this parameter is ignored.
    - For $N=2$, options are:
      - `0`: Periodic spline.
      - `1`: Starting derivative specified.
    - For $N=3$, options are:
      - `0`: Periodic spline.
      - `1`: Clamped spline.
      - `2`: Natural spline.
    - For $N \geq 4$, options are:
      - `0`: Periodic spline.
      - `1`: All starting derivatives specified.

  - **`boundary_condition`**: Not strictly enforced; only the first $N-1$ values need to match the required conditions.

The same parameter rules apply to the `Curve` fitting class.

---

For more details, refer to the `doxygen`-generated documentation or the project report in `doc/report.pdf`.
