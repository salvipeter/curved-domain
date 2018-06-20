#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <stack>
#include <vector>

#include <geometry.hh>
#include <harmonic.h>

using namespace Geometry;

const size_t LEVELS = 9;
const double DENSITY = 0.05;
const double EPSILON = 1.0e-5;

void floodFill(std::vector<bool> &grid, size_t n, size_t x, size_t y) {
  using Coord = std::pair<size_t, size_t>;
  std::stack<Coord> ps;
  ps.push({x, y});
  do {
    size_t x = ps.top().first, y = ps.top().second; // auto [x, y] = ps.top();
    ps.pop();
    if (grid[y*n+x])
      continue;
    grid[y*n+x] = true;
    if (x > 0)     ps.push({x - 1, y});
    if (x < n - 1) ps.push({x + 1, y});
    if (y > 0)     ps.push({x,     y - 1});
    if (y < n - 1) ps.push({x,     y + 1});
  } while (!ps.empty());
}

// Generates a mesh using a discretization of the domain (using a bitmap of size 2^size).
// Assumes that domain is in [-1,1]x[-1,1].
TriMesh regularMesh(const Point2DVector &domain, size_t size) {
  size_t n = (size_t)std::pow(2, size);
  std::vector<bool> grid(n * n, false);
  Point2D offset(-1.05, -1.05);
  double scaling = n / 2.1;

  // Init
  Point2D p0 = (domain.back() - offset) * scaling;
  int x0 = (int)p0[0], y0 = (int)p0[1];
  for (const auto &p : domain) {
    Point2D p1 = (p - offset) * scaling;
    int x1 = (int)p1[0], y1 = (int)p1[1];
    int dx = std::abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
    int dy = std::abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
    int err = (dx > dy ? dx : -dy) / 2, e2;
    while (true) {
      grid[y0*n+x0] = true;
      if (x0 == x1 && y0 == y1)
        break;
      e2 = err;
      if (e2 > -dx) { err -= dy; x0 += sx; }
      if (e2 <  dy) { err += dx; y0 += sy; }
    }
  }
  floodFill(grid, n, 0, 0);

  // Build mesh
  std::vector<size_t> row(n);
  PointVector pv;
  TriMesh result;
  size_t index = 0;

  for (size_t j = 1; j < n; ++j) {
    for (size_t i = 0; i < n - 1; ++i) {
      if (grid[j*n+i])
        continue;

      if (grid[(j-1)*n+i+1]) {
        // no NE
        if (!grid[(j-1)*n+i] && !grid[j*n+i+1])
          result.addTriangle(index, row[i], index + 1);   // N & E
      } else {
        // NE
        if (!grid[(j-1)*n+i])
          result.addTriangle(index, row[i], row[i+1]);    // N & NE
        if (!grid[j*n+i+1])
          result.addTriangle(index, row[i+1], index + 1); // E & NE
      }

      pv.emplace_back((double)i / scaling + offset[0], (double)j / scaling + offset[1], 0.0);
      row[i] = index++;
    }
  }

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

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " domain.dom" << std::endl;
    return 1;
  }

  // Read the domain from the given file
  // Format:
  // 4                     ; # of vertices (including concave corners)
  // 0 0 1 0               ; 0 for convex, 1 for concave
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
      } while(concave[k]);
      if (pv.size() == 6)
        harmonic_add_line(map, &pv[0], &pv[1]);
      else
        harmonic_add_curve(map, &pv[0], pv.size() / 3);
    }
    harmonic_solve(map, EPSILON, false);
    parameters.push_back(map);
  }

  // Compute the l,s,h values
  TriMesh mesh = regularMesh(points, LEVELS);
  std::vector<PointVector> values[3];
  for (size_t lhs = 0; lhs < 3; ++lhs)
    values[lhs].reserve(mesh.points().size());
  for (const auto &p : mesh.points()) {
    DoubleVector bc(parameters.size(), 0.0);
    for (size_t i = 0; i < n; ++i)
      harmonic_eval(parameters[i], p.data(), &bc[i]);
    PointVector l, h, s;
    for (size_t i = 0; i < parameters.size(); ++i) {
      Point2D sd = barycentricSD(bc, i);
      l.emplace_back(p[0], p[1], bc[i]);
      h.emplace_back(p[0], p[1], sd[1]);
      s.emplace_back(p[0], p[1], sd[0]);
    }
    values[0].push_back(l);
    values[1].push_back(h);
    values[2].push_back(s);
  }

  // Write the PS output
  {
    auto scale = [](Point2D p) { return (p + Point2D(1,1)) * 250 + Point2D(50,50); };
    std::ofstream f("curved-domain.eps");
    for (size_t lhs = 0; lhs < 3; ++lhs)
      for (size_t i = 0; i < parameters.size(); ++i) {
        f << "1 0 0 setrgbcolor\n"
          << "newpath\n";
        auto q = scale(points.back());
        f << q[0] << ' ' << q[1] << " moveto\n";
        for (const auto &p : points) {
          auto q = scale(p);
          f << q[0] << ' ' << q[1] << " lineto\n";
        }
        f << "stroke\n"
          << "0 0 0 setrgbcolor" << std::endl;
        writeContours(f, scale, mesh, values[lhs][i], DENSITY);
        f << "showpage" << std::endl;
      }
  }

  // Deallocate memory
  for (auto p : parameters)
    harmonic_free(p);

  return 0;
}
