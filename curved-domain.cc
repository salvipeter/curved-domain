#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <stack>
#include <vector>

#include <geometry.hh>
#include <harmonic.h>

using namespace Geometry;

const size_t LEVELS = 9;
const double EPSILON = 1.0e-5;
const size_t RESOLUTION = 50;   // Bezier curve resolution

TriMesh regularMesh(const Point2DVector &domain, size_t size) {
  size_t n = (size_t)std::pow(2, size);
  Point2D offset(-1.05, -1.05);
  double scaling = n / 2.1;

  PointVector pv;
  TriMesh result;

  for (size_t j = 0; j < n - 1; ++j) {
    for (size_t i = 0; i < n - 1; ++i) {
      result.addTriangle(i*n+j, i*n+j+1, (i+1)*n+j);
      result.addTriangle((i+1)*n+j, i*n+j+1, (i+1)*n+j+1);
    }
  }

  for (size_t j = 0; j < n; ++j)
    for (size_t i = 0; i < n; ++i)
      pv.emplace_back((double)i / scaling + offset[0], (double)j / scaling + offset[1], 0.0);

  result.setPoints(pv);
  return result;
}

void writeContours(std::ofstream &f, const std::function<Point2D(const Point2D &)> &scale,
                   const TriMesh &mesh, const PointVector &points, double density) {
  Point2DVector found;
  auto slice = [&points, density, &found](size_t i, size_t j) {
    double x = points[i][2], y = points[j][2];
    if (y > x) {
      std::swap(x, y);
      std::swap(i, j);
    }
    int q1 = static_cast<int>(std::floor(x / density));
    int q2 = static_cast<int>(std::floor(y / density));
    if (q1 - q2 == 1 && q1 != 0) {
      double alpha = (density * q1 - y) / (x - y);
      Point3D q = points[j] * (1.0 - alpha) + points[i] * alpha;
      found.emplace_back(q[0], q[1]);
    }
  };
  using Segment = std::pair<Point2D, Point2D>;
  std::vector<Segment> segments;
  for (const auto &tri : mesh.triangles()) {
    found.clear();
    slice(tri[0], tri[1]);
    slice(tri[0], tri[2]);
    slice(tri[1], tri[2]);
    if (found.size() == 2)
      segments.emplace_back(found[0], found[1]);
  }
  for (const auto &s : segments) {
    auto p = scale(s.first), q = scale(s.second);
    f << "newpath\n"
      << p[0] << ' ' << p[1] << " moveto\n"
      << q[0] << ' ' << q[1] << " lineto\n"
      << "stroke" << std::endl;
  }
}

// Rescale to [-1,-1]x[1,1].
void rescaleDomain(Point2DVector &domain) {
  double minx = 0.0, miny = 0.0, maxx = 0.0, maxy = 0.0;
  for (const auto &v : domain) {
    minx = std::min(minx, v[0]); miny = std::min(miny, v[1]);
    maxx = std::max(maxx, v[0]); maxy = std::max(maxy, v[1]);
  }
  double width = std::max(maxx - minx, maxy - miny);
  Point2D topleft(-1.0, -1.0), minp(minx, miny);
  for (auto &v : domain)
    v = topleft + (v - minp) * 2.0 / width;
}

// Returns the (s,d) system for side i, given the barycentric coordinates bc.
Point2D barycentricSD(const DoubleVector &bc, size_t i) {
  size_t n = bc.size(), im = (i + n - 1) % n;
  double him = bc[im], hi = bc[i];
  double d = 1.0 - him - hi;
  double s = him + hi;
  if (std::abs(s) > EPSILON)
    s = hi / s;
  return Point2D(s, d);
}

void bernsteinAll(size_t n, double u, DoubleVector &coeff) {
  coeff.clear(); coeff.reserve(n + 1);
  coeff.push_back(1.0);
  double u1 = 1.0 - u;
  for (size_t j = 1; j <= n; ++j) {
    double saved = 0.0;
    for (size_t k = 0; k < j; ++k) {
      double tmp = coeff[k];
      coeff[k] = saved + tmp * u1;
      saved = tmp * u;
    }
    coeff.push_back(saved);
  }
}

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " domain.dom [density]" << std::endl;
    std::cerr << "Density defaults to 0.1." << std::endl;
    return 1;
  }

  double density = 0.1;
  if (argc > 2)
    density = std::stod(argv[2]);

  // Read the domain from the given file
  // Format:
  // 4                     ; # of vertices (including concave corners)
  // 0 0 1 0               ; 0 for convex, 1 for concave [i.e., internal Bezier control point]
  // 0 0                   ; 1st vertex
  // 4 1                   ; 2nd vertex
  // 3 0                   ; 3rd vertex
  // 4 -1                  ; 4th vertex
  size_t n;
  Point2DVector points;
  std::vector<bool> concave;
  {
    std::ifstream f(argv[1]);
    f >> n;
    concave.resize(n);
    for (size_t i = 0; i < n; ++i) {
      int b;
      f >> b;
      concave[i] = b != 0;
    }
    double x, y;
    for (size_t i = 0; i < n; ++i) {
      f >> x >> y;
      points.emplace_back(x, y);
    }
    if (!f.good()) {
      std::cerr << "Cannot read file: " << argv[1] << std::endl;
      return 2;
    }
  }

  rescaleDomain(points);

  // Create parameterization
  std::vector<HarmonicMap *> parameters;
  DoubleVector min = { -1, -1 }, max = { 1, 1 };
  for (size_t i = 0; i < n; ++i) {
    if (concave[i])
      continue;
    auto map = harmonic_create(&min[0], &max[0], LEVELS);
    for (size_t j = 0; j < n; ++j) {
      if (concave[j])
        continue;
      DoubleVector pv;
      pv.push_back(points[j][0]); pv.push_back(points[j][1]); pv.push_back(j == i ? 1.0 : 0.0);
      size_t k = j;
      do {
        k = (k + 1) % n;
        pv.push_back(points[k][0]); pv.push_back(points[k][1]); pv.push_back(k == i ? 1.0 : 0.0);
      } while (concave[k]);
      if (pv.size() == 6)
        harmonic_add_line(map, &pv[0], &pv[3]);
      else {
        size_t m = pv.size() / 3;
        double from = pv[2], to = pv[m*3-1];
        if (from != to) {
          // Distribute the weights evenly
          for (size_t k = 1; k < m - 1; ++k) {
            double alpha = (double)k / (m - 1);
            pv[k*3+2] = std::max(from * (1.0 - alpha) + to * alpha, 0.0);
          }
        }
        harmonic_add_curve(map, &pv[0], m, RESOLUTION);
      }
    }
    harmonic_solve(map, EPSILON, false);
    parameters.push_back(map);
    // PPM output
    std::stringstream ppm;
    ppm << "harmonic" << i << ".ppm";
    // harmonic_write_ppm(map, ppm.str().c_str());
  }
  size_t m = parameters.size();

  // Compute the l,s,h values
  TriMesh mesh = regularMesh(points, LEVELS);
  std::vector<PointVector> values[3];
  for (size_t lhs = 0; lhs < 3; ++lhs) {
    values[lhs].resize(m);
    for (size_t i = 0; i < m; ++i)
      values[lhs][i].reserve(mesh.points().size());
  }
  for (const auto &p : mesh.points()) {
    DoubleVector bc(m, 0.0);
    for (size_t i = 0; i < m; ++i)
      harmonic_eval(parameters[i], p.data(), &bc[i]);
    for (size_t i = 0; i < m; ++i) {
      Point2D sd = barycentricSD(bc, i);
      values[0][i].emplace_back(p[0], p[1], bc[i]);
      values[1][i].emplace_back(p[0], p[1], sd[1]);
      values[2][i].emplace_back(p[0], p[1], sd[0]);
    }
  }

  // Write the PS output
  {
    auto scale = [](Point2D p) { return (p + Point2D(1,1)) * 250 + Point2D(50,50); };
    std::ofstream f("curved-domain.eps");
    for (size_t lhs = 0; lhs < 3; ++lhs)
      for (size_t i = 0; i < m; ++i) {
        // Draw the domain in red
        f << "1 0 0 setrgbcolor\n"
          << "newpath\n";
        bool first = true;
        for (size_t j = 0; j < n; ++j) {
          if (concave[j])
            continue;
          if (first) {
            first = false;
            auto q = scale(points[j]);
            f << q[0] << ' ' << q[1] << " moveto\n";
          }
          Point2DVector pv;
          pv.push_back(scale(points[j]));
          size_t k = j;
          do {
            k = (k + 1) % n;
            pv.push_back(scale(points[k]));            
          } while (concave[k]);
          if (pv.size() == 2) {
            // Line
            auto &q = pv.back();
            f << q[0] << ' ' << q[1] << " lineto\n";
          } else {
            // Bezier curve
            DoubleVector coeff;
            size_t d = pv.size() - 1;
            for (size_t k = 1; k <= RESOLUTION; ++k) {
              double u = (double)k / RESOLUTION;
              bernsteinAll(d, u, coeff);
              Point2D q(0.0, 0.0);
              for (size_t l = 0; l <= d; ++l)
                q += pv[l] * coeff[l];
              f << q[0] << ' ' << q[1] << " lineto\n";
            }
          }
        }
        f << "stroke\n";
        // Draw the contours in black
        f << "0 0 0 setrgbcolor" << std::endl;
        writeContours(f, scale, mesh, values[lhs][i], density);
        f << "showpage" << std::endl;
      }
  }

  // Deallocate memory
  for (auto p : parameters)
    harmonic_free(p);

  return 0;
}
