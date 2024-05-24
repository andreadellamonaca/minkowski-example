#include <iostream>
#include <boost/python.hpp>
#include <boost/polygon/polygon.hpp>
#include <vector>
#include <utility>
#include <chrono>

typedef boost::polygon::point_data<int> point;
typedef boost::polygon::polygon_set_data<int> polygon_set;
typedef boost::polygon::polygon_with_holes_data<int> polygon;
typedef std::pair<point, point> edge;
using namespace boost::polygon::operators;
using namespace boost::python;

void convolve_two_segments(std::vector<point>& figure, const edge& a, const edge& b) {
  using namespace boost::polygon;
  figure.clear();
  figure.push_back(point(a.first));
  figure.push_back(point(a.first));
  figure.push_back(point(a.second));
  figure.push_back(point(a.second));
  convolve(figure[0], b.second);
  convolve(figure[1], b.first);
  convolve(figure[2], b.first);
  convolve(figure[3], b.second);
}

template <typename itrT1, typename itrT2>
void convolve_two_point_sequences(polygon_set& result, itrT1 ab, itrT1 ae, itrT2 bb, itrT2 be) {
  using namespace boost::polygon;
  if(ab == ae || bb == be)
    return;
  point first_a = *ab;
  point prev_a = *ab;
  std::vector<point> vec;
  polygon poly;
  ++ab;
  for( ; ab != ae; ++ab) {
    point first_b = *bb;
    point prev_b = *bb;
    itrT2 tmpb = bb;
    ++tmpb;
    for( ; tmpb != be; ++tmpb) {
      convolve_two_segments(vec, std::make_pair(prev_b, *tmpb), std::make_pair(prev_a, *ab));
      set_points(poly, vec.begin(), vec.end());
      result.insert(poly);
      prev_b = *tmpb;
    }
    prev_a = *ab;
  }
}

template <typename itrT>
void convolve_point_sequence_with_polygons(polygon_set& result, itrT b, itrT e, const std::vector<polygon>& polygons) {
  using namespace boost::polygon;
  for(std::size_t i = 0; i < polygons.size(); ++i) {
    convolve_two_point_sequences(result, b, e, begin_points(polygons[i]), end_points(polygons[i]));
    for(polygon_with_holes_traits<polygon>::iterator_holes_type itrh = begin_holes(polygons[i]);
        itrh != end_holes(polygons[i]); ++itrh) {
      convolve_two_point_sequences(result, b, e, begin_points(*itrh), end_points(*itrh));
    }
  }
}

void convolve_two_polygon_sets(polygon_set& result, const polygon_set& a, const polygon_set& b) {
  using namespace boost::polygon;
  result.clear();
  std::vector<polygon> a_polygons;
  std::vector<polygon> b_polygons;
  a.get(a_polygons);
  b.get(b_polygons);
  for(std::size_t ai = 0; ai < a_polygons.size(); ++ai) {
    convolve_point_sequence_with_polygons(result, begin_points(a_polygons[ai]), 
                                          end_points(a_polygons[ai]), b_polygons);
    for(polygon_with_holes_traits<polygon>::iterator_holes_type itrh = begin_holes(a_polygons[ai]);
        itrh != end_holes(a_polygons[ai]); ++itrh) {
      convolve_point_sequence_with_polygons(result, begin_points(*itrh), 
                                            end_points(*itrh), b_polygons);
    }
    for(std::size_t bi = 0; bi < b_polygons.size(); ++bi) {
      polygon tmp_poly = a_polygons[ai];
      result.insert(convolve(tmp_poly, *(begin_points(b_polygons[bi]))));
      tmp_poly = b_polygons[bi];
      result.insert(convolve(tmp_poly, *(begin_points(a_polygons[ai]))));
    }
  }
}

namespace boost { namespace polygon{

  template <typename T>
  std::ostream& operator<<(std::ostream& o, const polygon_data<T>& poly) {
    o << "Polygon { ";
    for(typename polygon_data<T>::iterator_type itr = poly.begin(); 
        itr != poly.end(); ++itr) {
      if(itr != poly.begin()) o << ", ";
      o << (*itr).get(HORIZONTAL) << " " << (*itr).get(VERTICAL);
    } 
    o << " } ";
    return o;
  } 

  template <typename T>
  std::ostream& operator<<(std::ostream& o, const polygon_with_holes_data<T>& poly) {
    o << "Polygon With Holes { ";
    for(typename polygon_with_holes_data<T>::iterator_type itr = poly.begin(); 
        itr != poly.end(); ++itr) {
      if(itr != poly.begin()) o << ", ";
      o << (*itr).get(HORIZONTAL) << " " << (*itr).get(VERTICAL);
    } o << " { ";
    for(typename polygon_with_holes_data<T>::iterator_holes_type itr = poly.begin_holes();
        itr != poly.end_holes(); ++itr) {
      o << (*itr);
    }
    o << " } } ";
    return o;
  }
}}

// Function to parse polygon set and extract points
boost::python::list convertPolygon(const polygon_set& polySet) {
    boost::python::list resultList;
    std::vector<polygon> a_polygons;
    polySet.get(a_polygons);

    // Iterate over each polygon in the set
    for(std::size_t ai = 0; ai < a_polygons.size(); ++ai) {
    	// Extract vertices of the current polygon
        for (auto vertexIt = begin_points(a_polygons[ai]); vertexIt != end_points(a_polygons[ai]); ++vertexIt) {
            // Add coordinates to the result
            resultList.append(boost::python::make_tuple(vertexIt->x(), vertexIt->y()));
        }
    }

    return resultList;
}

polygon convertPythonList(const boost::python::list& pointsList) {
    std::vector<point> pts;
    polygon poly;

    for (int i = 0; i < boost::python::len(pointsList); ++i) {
        boost::python::tuple pointTuple = boost::python::extract<boost::python::tuple>(pointsList[i]);
        int x = boost::python::extract<int>(pointTuple[0]);
        int y = boost::python::extract<int>(pointTuple[1]);
        pts.push_back(point(x, y));
    }
    boost::polygon::set_points(poly, pts.begin(), pts.end());

    return poly;
}


boost::python::list minkowski_sum(boost::python::list &polyA, boost::python::list &polyB) {
  using std::chrono::high_resolution_clock;
  using std::chrono::duration_cast;
  using std::chrono::duration;
  using std::chrono::milliseconds;
  polygon_set a, b, c;
  
  std::vector<polygon> polys;
  std::vector<point> pts;
  auto t1 = high_resolution_clock::now();
  polygon poly = convertPythonList(polyA);
  a+=poly;
  auto t2 = high_resolution_clock::now();

  /* Getting number of milliseconds as an integer. */
  auto ms_int = duration_cast<milliseconds>(t2 - t1);
  std::cout << ms_int.count() << "ms\n";

  
  auto t3 = high_resolution_clock::now();
  polygon poly2 = convertPythonList(polyB);
  b+=poly2;
  auto t4 = high_resolution_clock::now();
  auto ms_int1 = duration_cast<milliseconds>(t4 - t3);
  std::cout << ms_int1.count() << "ms\n";

  polys.clear();
  auto t5 = high_resolution_clock::now();
  convolve_two_polygon_sets(c, a, b);
  auto t6 = high_resolution_clock::now();
  auto ms_int2 = duration_cast<milliseconds>(t6 - t5);
  std::cout << ms_int2.count() << "ms\n";

  auto t51 = high_resolution_clock::now();
  c.get(polys);
  auto t61 = high_resolution_clock::now();
  auto ms_int21 = duration_cast<milliseconds>(t61 - t51);
  std::cout << ms_int21.count() << "ms\n";

  auto t7 = high_resolution_clock::now();
  boost::python::list points = convertPolygon(c);
  auto t8 = high_resolution_clock::now();
  auto ms_int3 = duration_cast<milliseconds>(t8 - t7);
  std::cout << ms_int3.count() << "ms\n";
  
  return points;
}



BOOST_PYTHON_MODULE(minkowski) {
    def("minkowski_sum", minkowski_sum);
}


