/*
  GEO1015.2021
  hw03 
  --
  Author:           Daniël Dobson
  Student number:   5152739
  Author:           Katrin Meschin
  Student number:   5163889
  Author:           Chao Gao
  Student number:   5474493

*/

#include "GroundFilter.h"

// -- LAS reading and writing
#include <lasreader.hpp>
#include <laswriter.hpp>

// -- CGAL delaunay triangulation
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>

// -- CGAL kd-tree
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>

// -- CGAL objects
#include <CGAL/Vector_3.h>
#include <CGAL/Triangle_3.h>

// -- Other
#include <cmath>
typedef Kernel::Triangle_3 Triangle;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Search_traits_3<Kernel> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
typedef Neighbor_search::Tree Tree;


std::vector<std::vector<double>> make_bbox(const std::vector<Point>& pointcloud) {
    /*
    Function that takes a vector of points and extracts its bounding box.

    Input:
        pointcloud:     input point cloud (an Nx3 numpy array)
    Output:
        vector of vectors containing upperleft {x_min, y_min} and lowerright {x_max, y_max} corner.
     */
    auto xExtremes = std::minmax_element(pointcloud.begin(), pointcloud.end(),
                                         [](const Point& lhs, const Point& rhs) {
                                             return lhs.x() < rhs.x();
                                         });

    auto yExtremes = std::minmax_element(pointcloud.begin(), pointcloud.end(),
                                         [](const Point& lhs, const Point& rhs) {
                                             return lhs.y() < rhs.y();
                                         });

    std::vector<std::vector<double>> bbox;
    std::vector<double> upperLeft = {xExtremes.first->x(), yExtremes.first->y()};
    std::vector<double> lowerRight = {xExtremes.second->x(), yExtremes.second->y()};
    bbox = {upperLeft, lowerRight};

    return bbox;
}

std::vector<int> make_cells(const std::vector<std::vector<double>>& bbox, const double& resolution) {
    /*
    Function that takes a bounding box and desired resolution, to return the number of rows and cells for
    construction of a grid.

    Input:
        bbox:       vector containing all extremes of data's extent.
        resolution: of type double
    Output:
        vector of integers {rows, cols}
    */
    int CELLROWS = std::ceil((bbox[1][1] - bbox[0][1]) / resolution);
    int CELLCOLS = std::ceil((bbox[1][0] - bbox[0][0]) / resolution);
    std::vector<int> CELLROWSCOLS = {CELLROWS, CELLCOLS};
    return CELLROWSCOLS;
}

template < typename T>
std::pair<bool, int > findInVector(const std::vector<T>  & vecOfElements, const T  & element)
{
    std::pair<bool, int > result;
    // Find given element in vector
    auto it = std::find(vecOfElements.begin(), vecOfElements.end(), element);
    if (it != vecOfElements.end())
    {
        result.second = distance(vecOfElements.begin(), it);
        result.first = true;
    }
    else
    {
        result.first = false;
        result.second = -1;
    }
    return result;
}


void groundfilter_tin(const std::vector<Point>& pointcloud, const json& jparams) {
  /*
    Function that performs ground filtering using TIN refinement and writes the result to a new LAS file.

    Inputs:
      pointcloud: input point cloud (an Nx3 numpy array),
      jparams: a dictionary jparams with all the parameters that are to be used in this function:
        - resolution:    resolution (cellsize) for the initial grid that is computed as part of the ground filtering algorithm,
        - distance:      distance threshold used in the ground filtering algorithm,
        - angle:         angle threshold used in the ground filtering algorithm in degrees,
        - output_las:    path to output .las file that contains your ground classification,
  */
  typedef CGAL::Projection_traits_xy_3<Kernel>  Gt;
  typedef CGAL::Delaunay_triangulation_2<Gt> DT;

  double resolution = jparams["resolution"];
  double distance = jparams["distance"];
  double angle = jparams["angle"];
  std::string output_las = jparams["output_las"];

  // Initialize bbox and number of cells blocks corresponding to resolution
  const std::vector<std::vector<double>> bbox = make_bbox(pointcloud);
  const std::vector<int> CELLROWSCOLS = make_cells(bbox, resolution);
  const int CELLROWS = CELLROWSCOLS[0];
  const int CELLCOLS = CELLROWSCOLS[1];

  // The offset is needed for determining the correct spacing of the grid cells
  const double offsetX = bbox[0][0];
  const double offsetY = bbox[0][1];


  int stabilityCount = 0; //
  while (stabilityCount < 1) { // this is an optional method for stability testing. Increase number larger than 1 to run more than once.
      // Initialize virtual grid to store initial ground points
      std::vector<std::vector<Point>> vGrid(CELLROWS, std::vector<Point>(CELLCOLS));

      for (int i = 0; i < CELLROWS; i++) {
          for (int j = 0; j < CELLCOLS; j++) {
              vGrid[i][j] = Point(0.,0.,0.);

          }
      }

      // Initialise vector to store class labels, codes used:
      // ‘2’ for ground points
      // ‘1’ for all the other points
      std::vector<int> class_labels;

      // For the construction of the rudimentary initial Delaunay TIN, the locally lowest elevation points are needed.
      // To achieve this, we loop over all points p in pointcloud, determine in which cellblock it would fall,
      // and only save the lowest local z-value for each cellblock.

      //std::vector<Point>::const_iterator itPc; // variable needed for finding index
      for (auto p: pointcloud) {

          // Determination to which cellblock point p belongs
          int cellX, cellY;
          cellX = std::floor((p.y() - offsetY) / resolution);
          cellY = std::floor((p.x() - offsetX) / resolution);


          // If grid cell is empty, we assign point p to it and the corresponding point p is marked as ground point.
          if ((vGrid[cellX][cellY][0] == 0.) && (vGrid[cellX][cellY][1] == 0.) && (vGrid[cellX][cellY][2] == 0.)){
              vGrid[cellX][cellY] = p;
              class_labels.push_back(2);
          }

              // If p has a lower elevation than previous lowest local elevation, grid cell is overwritten with p.
              // The previous point's class is reverted to non-ground point, and current point is marked as ground point.
          else if (p[2] < vGrid[cellX][cellY][2]) {

              std::pair<bool, int> result = findInVector<Point>(pointcloud, vGrid[cellX][cellY]);
              if (result.first)
                  int nothing;
              else {
                  // Useful for debugging
                  std::cout << "Element Not Found" << std::endl;
                  std::cout << "p:          " << p << "\n";
                  std::cout << "vGrid[cellX][cellY]:    " << vGrid[cellX][cellY] << "\n";
              }
              class_labels[result.second] = 1;
              vGrid[cellX][cellY] = p;
              class_labels.push_back(2);
          }

              // If neither, the point is skipped and marked as non-ground point.
          else {
              class_labels.push_back(1);
          }

      }


      // Initialise the Delaunay Triangulation (DT) object
      DT dt;

      // Insertion of initial ground points from vGrid into DT
      for (int i = 0; i < CELLROWS; i++) {
          for (int j = 0; j < CELLCOLS; j++) {
              dt.insert(vGrid[i][j]);
          }
      }

      // Computation of geometric properties for all points p in pointcloud:
      // 1. distance:   Orthogonal point-plane distance of point p and corresponding triangle t from Dt.
      // 2. alpha:      largest angle of angles between point p and the vectors that connect it to corresponding triangle t
      // vertices.
      int countComp = 0;
      for (auto p: pointcloud) {
          if (class_labels[countComp] != 2) {
              DT::Face_handle triangle = dt.locate(p);

              // The 3 vertices of the triangle:
              DT::Vertex_handle v0 = triangle->vertex(0);
              DT::Vertex_handle v1 = triangle->vertex(1);
              DT::Vertex_handle v2 = triangle->vertex(2);

              // Put points in appropriate triangle/plane object
              Triangle located;
              located = Triangle(v0->point(), v1->point(), v2->point());

              // Orthogonal distance calculation
              double dist = std::sqrt(CGAL::squared_distance(p, located));

              if (dist < distance) {
                  // Max angle calculation
                  std::vector<double> betas;
                  double beta1, beta2, beta3, alphamax;
                  //alphamax = 0.;
                  Point p0;
                  Vector pv0, pv1, pv2, p0v0, p0v1, p0v2;



                  pv0 = Vector(p, v0->point());
                  pv1 = Vector(p, v1->point()), pv2 = Vector(p, v2->point());
                  p0 = Point(p[0], p[1], p[2] - dist);
                  p0v0 = Vector(p0, v0->point());
                  p0v1 = Vector(p0, v1->point()), p0v2 = Vector(p0, v2->point());
                  beta1 = std::acos((pv0 * p0v0) / ((std::sqrt(pv0.squared_length()) * std::sqrt(p0v0.squared_length()))));
                  betas.push_back(beta1);
                  beta2 = std::acos((pv1 * p0v1) / ((std::sqrt(pv1.squared_length()) * std::sqrt(p0v1.squared_length()))));
                  betas.push_back(beta2);
                  beta3 = std::acos((pv2 * p0v2) / ((std::sqrt(pv2.squared_length()) * std::sqrt(p0v2.squared_length()))));
                  betas.push_back(beta3);
                  auto it = std::max_element(betas.begin(), betas.end());
                  alphamax = it[0] * 180 / M_PI;


                  if (alphamax < angle) {
                      dt.insert(p);
                      class_labels[countComp] = 2;
                  }
              }
          }
          countComp++;
      }

      // Print number of ground points to console
      int groundPoints = 0, otherPoints = 0;
      for (auto label: class_labels) { if (label == 2) { groundPoints++; } else { otherPoints++; }}
      std::cout << "Number of classified ground points:         " << groundPoints << "\n";
      std::cout << "Number of other (non-ground) points:        " << otherPoints << "\n";
      std::cout << "Sum:                                        " << groundPoints + otherPoints << "\n";
      std::cout << "Original size pointcloud:                   " << pointcloud.size() << "\n";

      // Write results to new LAS file
      write_lasfile(jparams["output_las"], pointcloud, class_labels);

      if (groundPoints == 0) {
          std::cout << "ERROR: No ground points detected.\n Run again, and/or change parameters.\n";
      }


      stabilityCount++;
  }

}

void groundfilter_csf(const std::vector<Point>& pointcloud, const json& jparams) {
  /*
  Function that performs ground filtering using CSF and writes the result to a new LAS file.

  Inputs:
    pointcloud: input point cloud (an Nx3 numpy array),
    jparams: a dictionary with all the parameters that are to be used in this function:
      - resolution:     resolution of the cloth grid,
      - epsilon_zmax:   tolerance to stop the iterations,
      - epsilon_ground: threshold used to classify ground points,
      - output_las:     path to output .las file that contains your ground classification
  */

  double resolution = jparams["resolution"];
  double epsilon_zmax = jparams["epsilon_zmax"];
  double epsilon_ground = jparams["epsilon_ground"];
  std::string output_las = jparams["output_las"];
  double displacement = 4;
  double scaling = 1.;

  // Create S inverted 3D and S 2D (for later use in kd-tree query).
  std::vector<Point> S3Dinverted;
  std::vector<Point> S2D;
  for (auto p : pointcloud) {
      Point p3Dinv = Point(p.x(), p.y(), p.z() * - 1);
      Point p2D = Point(p.x(), p.y(), 0.);
      S3Dinverted.push_back(p3Dinv);
      S2D.push_back(p2D);
  }

  // Initialise the cloth C at an elevation z0 higher than the highest elevation
  double maxZ;
  for (auto p : S3Dinverted) {if (p.z() > maxZ) maxZ = p.z();}
  double z0 = maxZ + (4 * displacement);

  // Initialize bbox and number of cells blocks for grid (cloth) corresponding to resolution
  const std::vector<std::vector<double>> bbox = make_bbox(pointcloud);
  const std::vector<int> CELLROWSCOLS = make_cells(bbox, resolution);
  const int CELLROWS = CELLROWSCOLS[0];
  const int CELLCOLS = CELLROWSCOLS[1];

  const double offsetX = bbox[0][0];
  const double offsetY = bbox[0][1];

  // Initialize cloth as regularly distributed grid of particles
  std::vector<std::vector<Point>> cloth(CELLROWS,std::vector<Point>(CELLCOLS));

  // Assign cloth points appropriate Point values
  for (int i=0; i < CELLROWS; i++) {
      for (int j=0; j < CELLCOLS; j++) {
          double cellX, cellY, cellZ;
          cellX = offsetX + (j * resolution);
          cellY = offsetY + (i * resolution);
          cellZ = z0;
          cloth[i][j] = Point(cellX, cellY, cellZ);
      }
  }


  // Initialise vector to store z parameters
  std::vector<std::vector<std::vector<double>>> zParams(CELLROWS,std::vector<std::vector<double>>(CELLCOLS));

  // Calculate particle elevation parameters
  // Construct and query kd-tree:
  Tree tree(S2D.begin(), S2D.end());
  const unsigned int N = 1;
  std::vector<Point>::iterator itS2D;

  for (int i=0; i < CELLROWS; i++) {
      for (int j=0; j < CELLCOLS; j++) {

          double pZmin, pZprev, pZcur;
          std::vector<double> zParam;
          Point query_point = Point(cloth[i][j].x(), cloth[i][j].y(), 0.);
          Neighbor_search search_result(tree, query_point, N);
          for (auto res: search_result) {
              Point neighbour_point = res.first;
              itS2D = std::find(S2D.begin(), S2D.end(), neighbour_point);
              int idxS3D = std::distance(S2D.begin(), itS2D);
              pZmin = S3Dinverted[idxS3D].z();
              zParams[i][j].push_back(pZmin);
              pZcur = z0;
              pZprev = z0 + displacement;
              zParams[i][j].push_back(pZprev);
              zParams[i][j].push_back(pZcur);
              zParams[i][j].push_back(1.); // Insert movable variable
          }
      }
  }

  double deltaZ = displacement;
  int deltaZstopper = 0;
  while (deltaZ > epsilon_zmax) {

      // external forces, apply to all movable p
      for (int i = 0; i < CELLROWS; i++) {
          for (int j = 0; j < CELLCOLS; j++) {
              if (zParams[i][j][3] == 1.) {
                  double pzmin = zParams[i][j][0];
                  double pzprev = zParams[i][j][1];
                  double pzcur = zParams[i][j][2];
                  double tmp = pzcur;
                  pzcur = (pzcur - pzprev) + pzcur;
                  pzprev = tmp;
                  if (pzcur <= pzmin) {pzcur = pzmin; zParams[i][j][3] = 0.;}
                  zParams[i][j][1] = pzprev;
                  zParams[i][j][2] = pzcur;
              }
          }
      }

      // internal forces, process once each set neighbours e of adjacent particles
      for (int i = 0; i < CELLROWS; i++) {
          for (int j = 0; j < CELLCOLS; j++) {
              if (zParams[i][j][3] == 1.) {
                  if ((i > 0) && (i < (CELLROWS - 1)) && (j > 0) && (j < (CELLCOLS - 1))) {

                      // n1
                      if (zParams[i-1][j][3] == 0.){
                          if (zParams[i-1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                          }
                          else if (zParams[i-1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                          }
                      }
                      if (zParams[i-1][j][3] == 1.){
                          if (zParams[i-1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                              zParams[i-1][j][2] += scaling * displacement;
                          }
                          else if (zParams[i-1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                              zParams[i-1][j][2] -= scaling * displacement;
                          }
                      }
                      // n2
                      if (zParams[i+1][j][3] == 0.){
                          if (zParams[i+1][j][2] < cloth[i][j].z()){
                          cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                        }
                          else if (zParams[i+1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                          }
                      }
                      if (zParams[i+1][j][3] == 1.){
                          if (zParams[i+1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                              zParams[i+1][j][2] += scaling * displacement;
                          }
                          else if (zParams[i+1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                              zParams[i+1][j][2] -= scaling * displacement;
                          }
                      }
                      // n3
                      if (zParams[i][j-1][3] == 0.){
                          if (zParams[i][j-1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                          }
                          else if (zParams[i][j-1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                          }
                      }
                      if (zParams[i][j-1][3] == 1.){
                          if (zParams[i][j-1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                              zParams[i][j-1][2] += scaling * displacement;
                          }
                          else if (zParams[i][j-1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                              zParams[i][j-1][2] -= scaling * displacement;
                          }
                      }
                      //n4
                      if (zParams[i][j+1][3] == 0.){
                          if (zParams[i][j+1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                          }
                          else if (zParams[i][j+1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                          }
                      }
                      if (zParams[i][j+1][3] == 1.){
                          if (zParams[i][j+1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                              zParams[i][j+1][2] += scaling * displacement;
                          }
                          else if (zParams[i][j+1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                              zParams[i][j+1][2] -= scaling * displacement;
                          }
                      }

                  } // 5
                  else if ((i == 0) && (j > 0) && (j < (CELLCOLS - 1))) {

                      // n1
                      if (zParams[i][j-1][3] == 0.){
                          if (zParams[i][j-1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                          }
                          else if (zParams[i][j-1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                          }
                      }
                      if (zParams[i][j-1][3] == 1.){
                          if (zParams[i][j-1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                              zParams[i][j-1][2] += scaling * displacement;
                          }
                          else if (zParams[i][j-1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                              zParams[i][j-1][2] -= scaling * displacement;
                          }
                      }
                      // n2
                      if (zParams[i][j+1][3] == 0.){
                          if (zParams[i][j+1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                          }
                          else if (zParams[i][j+1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                          }
                      }
                      if (zParams[i][j+1][3] == 1.){
                          if (zParams[i][j+1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                              zParams[i][j+1][2] += scaling * displacement;
                          }
                          else if (zParams[i][j+1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                              zParams[i][j+1][2] -= scaling * displacement;
                          }
                      }
                      // n3
                      if (zParams[i+1][j][3] == 0.){
                          if (zParams[i+1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                          }
                          else if (zParams[i+1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                          }
                      }
                      if (zParams[i+1][j][3] == 1.){
                          if (zParams[i+1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                              zParams[i+1][j][2] += scaling * displacement;
                          }
                          else if (zParams[i+1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                              zParams[i+1][j][2] -= scaling * displacement;
                          }
                      }

                  } // 2
                  else if ((i > 0) && (i < (CELLROWS - 1)) && (j == 0)) {

                      // n1
                      if (zParams[i][j+1][3] == 0.){
                          if (zParams[i][j+1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                          }
                          else if (zParams[i][j+1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                          }
                      }
                      if (zParams[i][j+1][3] == 1.){
                          if (zParams[i][j+1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                              zParams[i][j+1][2] += scaling * displacement;
                          }
                          else if (zParams[i][j+1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                              zParams[i][j+1][2] -= scaling * displacement;
                          }
                      }
                      // n2
                      if (zParams[i-1][j][3] == 0.){
                          if (zParams[i-1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                          }
                          else if (zParams[i-1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                          }
                      }
                      if (zParams[i-1][j][3] == 1.){
                          if (zParams[i-1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                              zParams[i-1][j][2] += scaling * displacement;
                          }
                          else if (zParams[i-1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                              zParams[i-1][j][2] -= scaling * displacement;
                          }
                      }
                      // n3
                      if (zParams[i+1][j][3] == 0.){
                          if (zParams[i+1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                          }
                          else if (zParams[i+1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                          }
                      }
                      if (zParams[i+1][j][3] == 1.){
                          if (zParams[i+1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                              zParams[i+1][j][2] += scaling * displacement;
                          }
                          else if (zParams[i+1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                              zParams[i+1][j][2] -= scaling * displacement;
                          }
                      }

                  } // 4
                  else if ((i > 0) && (i < (CELLROWS - 1)) && (j == (CELLCOLS - 1))) {

                      // n1
                      if (zParams[i][j-1][3] == 0.){
                          if (zParams[i][j-1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                          }
                          else if (zParams[i][j-1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                          }
                      }
                      if (zParams[i][j-1][3] == 1.){
                          if (zParams[i][j-1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                              zParams[i][j-1][2] += scaling * displacement;
                          }
                          else if (zParams[i][j-1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                              zParams[i][j-1][2] -= scaling * displacement;
                          }
                      }
                      // n2
                      if (zParams[i-1][j][3] == 0.){
                          if (zParams[i-1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                          }
                          else if (zParams[i-1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                          }
                      }
                      if (zParams[i-1][j][3] == 1.){
                          if (zParams[i-1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                              zParams[i-1][j][2] += scaling * displacement;
                          }
                          else if (zParams[i-1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                              zParams[i-1][j][2] -= scaling * displacement;
                          }
                      }
                      //n3
                      if (zParams[i+1][j][3] == 0.){
                          if (zParams[i+1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                          }
                          else if (zParams[i+1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                          }
                      }
                      if (zParams[i+1][j][3] == 1.){
                          if (zParams[i+1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                              zParams[i+1][j][2] += scaling * displacement;
                          }
                          else if (zParams[i+1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                              zParams[i+1][j][2] -= scaling * displacement;
                          }
                      }

                  } // 6
                  else if ((i == (CELLROWS - 1)) && (j > 0) && (j < (CELLCOLS - 1))) {

                      // n1
                      if (zParams[i][j-1][3] == 0.){
                          if (zParams[i][j-1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                          }
                          else if (zParams[i][j-1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                          }
                      }
                      if (zParams[i][j-1][3] == 1.){
                          if (zParams[i][j-1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                              zParams[i][j-1][2] += scaling * displacement;
                          }
                          else if (zParams[i][j-1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                              zParams[i][j-1][2] -= scaling * displacement;
                          }
                      }
                      // n2
                      if (zParams[i][j+1][3] == 0.){
                          if (zParams[i][j+1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                          }
                          else if (zParams[i][j+1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                          }
                      }
                      if (zParams[i][j+1][3] == 1.){
                          if (zParams[i][j+1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                              zParams[i][j+1][2] += scaling * displacement;
                          }
                          else if (zParams[i][j+1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                              zParams[i][j+1][2] -= scaling * displacement;
                          }
                      }
                      // n3
                      if (zParams[i-1][j][3] == 0.){
                          if (zParams[i-1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                          }
                          else if (zParams[i-1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                          }
                      }
                      if (zParams[i-1][j][3] == 1.){
                          if (zParams[i-1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                              zParams[i-1][j][2] += scaling * displacement;
                          }
                          else if (zParams[i-1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                              zParams[i-1][j][2] -= scaling * displacement;
                          }
                      }
                  } // 8
                  else if ((i == 0) && (j == 0)) {

                      // n1
                      if (zParams[i][j+1][3] == 0.){
                          if (zParams[i][j+1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                          }
                          else if (zParams[i][j+1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                          }
                      }
                      if (zParams[i][j+1][3] == 1.){
                          if (zParams[i][j+1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                              zParams[i][j+1][2] += scaling * displacement;
                          }
                          else if (zParams[i][j+1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                              zParams[i][j+1][2] -= scaling * displacement;
                          }
                      }
                      // n2
                      if (zParams[i+1][j][3] == 0.){
                          if (zParams[i+1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                          }
                          else if (zParams[i+1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                          }
                      }
                      if (zParams[i+1][j][3] == 1.){
                          if (zParams[i+1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                              zParams[i+1][j][2] += scaling * displacement;
                          }
                          else if (zParams[i+1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                              zParams[i+1][j][2] -= scaling * displacement;
                          }
                      }

                  } // 1
                  else if ((i == 0) && (j == (CELLCOLS - 1))) {

                      // n1
                      if (zParams[i][j-1][3] == 0.){
                          if (zParams[i][j-1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                          }
                          else if (zParams[i][j-1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                          }
                      }
                      if (zParams[i][j-1][3] == 1.){
                          if (zParams[i][j-1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                              zParams[i][j-1][2] += scaling * displacement;
                          }
                          else if (zParams[i][j-1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                              zParams[i][j-1][2] -= scaling * displacement;
                          }
                      }
                      // n2
                      if (zParams[i+1][j][3] == 0.){
                          if (zParams[i+1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                          }
                          else if (zParams[i+1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                          }
                      }
                      if (zParams[i+1][j][3] == 1.){
                          if (zParams[i+1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                              zParams[i+1][j][2] += scaling * displacement;
                          }
                          else if (zParams[i+1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                              zParams[i+1][j][2] -= scaling * displacement;
                          }
                      }

                  } // 3
                  else if ((i == (CELLROWS - 1)) && (j == 0)) {

                      // n1
                      if (zParams[i-1][j][3] == 0.){
                          if (zParams[i-1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                          }
                          else if (zParams[i-1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                          }
                      }
                      if (zParams[i-1][j][3] == 1.){
                          if (zParams[i-1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                              zParams[i-1][j][2] += scaling * displacement;
                          }
                          else if (zParams[i-1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                              zParams[i-1][j][2] -= scaling * displacement;
                          }
                      }
                      // n2
                      if (zParams[i][j+1][3] == 0.){
                          if (zParams[i][j+1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                          }
                          else if (zParams[i][j+1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                          }
                      }
                      if (zParams[i][j+1][3] == 1.){
                          if (zParams[i][j+1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                              zParams[i][j+1][2] += scaling * displacement;
                          }
                          else if (zParams[i][j+1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                              zParams[i][j+1][2] -= scaling * displacement;
                          }
                      }

                  } // 7
                  else if ((i == (CELLROWS - 1)) && (j == (CELLCOLS - 1))) {

                      // n1
                      if (zParams[i][j-1][3] == 0.){
                          if (zParams[i][j-1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                          }
                          else if (zParams[i][j-1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                          }
                      }
                      if (zParams[i][j-1][3] == 1.){
                          if (zParams[i][j-1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                              zParams[i][j-1][2] += scaling * displacement;
                          }
                          else if (zParams[i][j-1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                              zParams[i][j-1][2] -= scaling * displacement;
                          }
                      }
                      // n2
                      if (zParams[i-1][j][3] == 0.){
                          if (zParams[i-1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                          }
                          else if (zParams[i-1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                          }
                      }
                      if (zParams[i-1][j][3] == 1.){
                          if (zParams[i-1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - scaling * displacement);
                              zParams[i-1][j][2] += scaling * displacement;
                          }
                          else if (zParams[i-1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + scaling * displacement);
                              zParams[i-1][j][2] -= scaling * displacement;
                          }
                      }

                  } // 9
              }
          }
      }
      // Calculate the max deltaZ
      for (int i=0; i < CELLROWS; i++) {
          for (int j=0; j < CELLCOLS; j++) {
              if ((zParams[i][j][2] - zParams[i][j][1]) > deltaZ) {deltaZ = zParams[i][j][2] - zParams[i][j][1];}
              if ((zParams[i][j][2] - zParams[i][j][1]) == deltaZ){deltaZstopper++;}
          }
      }
      if (deltaZstopper >= 500){break;}
  }

  std::vector<Point> clothVec;
    for (int i=0; i < CELLROWS; i++) {
        for (int j=0; j < CELLCOLS; j++) {
            clothVec.push_back(cloth[i][j]);
        }
    }

  // Initialise vector to store class labels, codes used:
  // ‘2’ for ground points
  // ‘1’ for all the other points
  std::vector<int> class_labels;


  Tree tree3d(clothVec.begin(), clothVec.end());
  for (auto p : S3Dinverted) {
      Point query_point = p;
      Neighbor_search search_result(tree3d, query_point, N);
      for (auto res : search_result) {
          Point neighbour_point = res.first;
          double distance = std::sqrt(res.second);
          if (distance < epsilon_ground) {class_labels.push_back(2);}
          else {class_labels.push_back(1);}
      }
  }

  // Print number of ground points to console
  int groundPoints=0, otherPoints=0;
  for (auto label : class_labels) {if (label == 2){groundPoints++;} else {otherPoints++;}}
  std::cout << "Number of classified ground points:         " << groundPoints << "\n";
  std::cout << "Number of other (non-ground) points:        " << otherPoints << "\n";
  std::cout << "Sum:                                        " << groundPoints + otherPoints << "\n";
  std::cout << "Original size pointcloud:                   " << pointcloud.size() << "\n";

  // Write the results to a new LAS file
  write_lasfile(jparams["output_las"], pointcloud, class_labels);

  if (groundPoints == 0) {
      std::cout << "WARNING: No ground points detected.\n Run again, and/or change parameters.\n";
  }
}



std::vector<Point> read_lasfile(const json& jparams) {
  /*
  Function to read points from a LAS file

  Inputs:
    jparams["filename"]:   the filename to read the LAS file to

  Returns:
    a std::vector<Point> with the points from the LAS file
  */
  std::string filename = jparams["filename"];
	LASreadOpener lasreadopener;
	lasreadopener.set_file_name(filename.c_str());
	LASreader* lasreader = lasreadopener.open();
	
	if (!lasreader){
		std::cerr << "cannot read las file: " << filename << "\n";
		exit(1);
	}

  //-- store each point in a CGAL Point_3 object
  //-- https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Point__3.html
	std::vector<Point> points;
	while (lasreader->read_point()) {
		points.push_back( 
			Point(
				lasreader->point.get_x(),
				lasreader->point.get_y(),
				lasreader->point.get_z()
			)
		);
	}
	lasreader->close();
	delete lasreader;

	return points;
}


void write_lasfile(const std::string filename, const std::vector<Point>& pointcloud, const std::vector<int>& class_labels) {
  /*
  Function to write a new LAS file with point labels (for the LAS classification field)

  Inputs:
    filename:   the filename to write the LAS file to
    pointcloud: input point cloud (a vector of Points),
    Labels:     Contains point labels. Should be a vector of ints of the same size as pointcloud (ie. one label for each point in the same order as pointcloud). Uses LAS classification codes, ie 2 = ground. 1 = unclassified.
  */
  LASwriteOpener laswriteopener;
  laswriteopener.set_file_name(filename.c_str());

  LASheader lasheader;
  lasheader.x_scale_factor = 0.01;
  lasheader.y_scale_factor = 0.01;
  lasheader.z_scale_factor = 0.01;
  lasheader.x_offset = 0.0;
  lasheader.y_offset = 0.0;
  lasheader.z_offset = 0.0;
  lasheader.point_data_format = 0;
  lasheader.point_data_record_length = 20;

  LASpoint laspoint;
  laspoint.init(&lasheader, lasheader.point_data_format, lasheader.point_data_record_length, 0);

  LASwriter* laswriter = laswriteopener.open(&lasheader);
  if (laswriter == 0)
  {
    std::cerr << "ERROR: could not open laswriter\n";
    exit(1);
  }

	if (pointcloud.size()!=class_labels.size()) {
		std::cerr << "ERROR: points has a different size than class_labels\n";
		exit(1);
	}

  for (size_t i=0; i<pointcloud.size(); ++i) {
		const Point& p = pointcloud[i];
		const int& label = class_labels[i];

    laspoint.set_x(p[0]);
    laspoint.set_y(p[1]);
    laspoint.set_z(p[2]);
		laspoint.set_classification(label);

    laswriter->write_point(&laspoint);
    laswriter->update_inventory(&laspoint);    
  } 

  laswriter->update_header(&lasheader, TRUE);
  laswriter->close();
  delete laswriter;
}