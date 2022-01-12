/*
  GEO1015.2021
  hw03 
  --
  Author:           Daniël Dobson
  Student number:   5152739

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
    int CELLROWS = std::ceil((bbox[1][1] - bbox[0][1]) / resolution);
    int CELLCOLS = std::ceil((bbox[1][0] - bbox[0][0]) / resolution);
    std::vector<int> CELLROWSCOLS = {CELLROWS, CELLCOLS};
    return CELLROWSCOLS;
}

void displace_pt(std::vector<double>& clothPt, std::vector<double> displacement) {
    std::transform(clothPt.begin(), clothPt.end(),
                   displacement.begin(), clothPt.begin(),
                   std::minus<double>());
}


void groundfilter_tin(const std::vector<Point>& pointcloud, const json& jparams) {
  /*
    !!! TO BE COMPLETED !!!
      
    Function that performs ground filtering using TIN refinement and writes the result to a new LAS file.

    !!! You are free to subdivide the functionality of this function into several functions !!!
      
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

  const double offsetX = bbox[0][0];
  const double offsetY = bbox[0][1];

  std::cout << "minX:     " << offsetX << "\n";
  std::cout << "minY:     " << offsetY << "\n";
  std::cout << "maxX:     " << bbox[1][0] << "\n";
  std::cout << "maxY:     " << bbox[1][1] << "\n";

  // Initialize virtual grid to store initial ground points
  std::vector<std::vector<Point>> vGrid(CELLROWS,std::vector<Point>(CELLCOLS));


  // Initialise vector to store class labels, codes used:
  // ‘2’ for ground points
  // ‘1’ for all the other points

  std::vector<int> class_labels;

  // The virtual grid will have one point per cellblock, which will contain lowest elevation value.

  std::vector<Point>::const_iterator itPc; // variable needed for finding index
  for (auto p : pointcloud) {
      int cellX, cellY;
      cellX = std::round((p.y() - offsetY) / resolution);
      cellY = std::round((p.x() - offsetX) / resolution);
      //std::cout << "x, y: " << cellX << ", " << cellY << std::endl;
      //std::cout << "what's inside this cell: " << vGrid[cellX][cellY] << std::endl;

      if ((vGrid[cellX][cellY][0] == 0. && vGrid[cellX][cellY][1] == 0. && vGrid[cellX][cellY][2] == 0.)) {
          vGrid[cellX][cellY] = p;
          class_labels.push_back(2);

      }
      else if (p[2] < vGrid[cellX][cellY][2]) {

          itPc = std::find(pointcloud.begin(), pointcloud.end(), vGrid[cellX][cellY]);

          int idxPc = std::distance(pointcloud.begin(), itPc);

          class_labels[idxPc] = 1;

          vGrid[cellX][cellY] = p;
          class_labels.push_back(2);
      }
      else {
          class_labels.push_back(1);
      }
  }


  // Initialise the Delaunay Triangulation (DT) object
  DT dt;

  // Insertion of initial points from vGrid into DT
  for (int i=0; i < CELLROWS; i++) {
      for (int j=0; j < CELLCOLS; j++) {
          dt.insert(vGrid[i][j]);
      }
  }

  // Computation of property distance
  int countComp = 0;
  for (auto p : pointcloud) {
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

          // Max angle calculation
          std::vector<double> betas;
          double beta1, beta2, beta3, betamax;
          Point p0;
          Vector pv0, pv1, pv2, p0v0, p0v1, p0v2;
          pv0 = Vector(p, v0->point()); pv1 = Vector(p, v1->point()), pv2 = Vector(p, v2->point());
          p0 = Point(p[0], p[1], p[2] - dist);
          p0v0 = Vector(p0, v0->point()); p0v1 = Vector(p0, v1->point()), p0v2 = Vector(p0, v2->point());
          beta1 = std::acos((pv0 * p0v0) / ((std::sqrt(pv0.squared_length()) * std::sqrt(p0v0.squared_length()))));
          betas.push_back(beta1);
          beta2 = std::acos((pv1 * p0v1) / ((std::sqrt(pv1.squared_length()) * std::sqrt(p0v1.squared_length()))));
          betas.push_back(beta2);
          beta3 = std::acos((pv2 * p0v2) / ((std::sqrt(pv2.squared_length()) * std::sqrt(p0v2.squared_length()))));
          betas.push_back(beta3);
          auto it = std::max_element(betas.begin(), betas.end());
          betamax = it[0]*180/M_PI;

          if ((dist < distance) && (betamax < angle)) {dt.insert(p); class_labels[countComp] = 2;}
      }
      countComp++;
  }

  // Write results to new LAS file
  write_lasfile(jparams["output_las"], pointcloud, class_labels);
}

void groundfilter_csf(const std::vector<Point>& pointcloud, const json& jparams) {
  /*
  !!! TO BE COMPLETED !!!
    
  Function that performs ground filtering using CSF and writes the result to a new LAS file.

  !!! You are free to subdivide the functionality of this function into several functions !!!
    
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

  //-- print the first 5 points in the pointcloud, which are CGAL Point_3
  //-- https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Point__3.html

  // Create S inverse 3D and S inverse 2D (for later use in kd-tree query)
  std::vector<Point> S3Dinverted;
  std::vector<Point> S2D;
  for (auto p : pointcloud) {
      Point p3Dinv = Point(p.x(), p.y(), p.z() * - 1);
      Point p2D = Point(p.x(), p.y(), 0.);
      S3Dinverted.push_back(p3Dinv);
      S2D.push_back(p2D);
  }

  //std::cout << "size pointcloud:        " << pointcloud.size() << std::endl;
  //std::cout << "size S3Dinverted: " << S3Dinverted.size() << std::endl;

  // Initialise the cloth C at an elevation z0 higher than the highest elevation
  double maxZ;
  double z0 = 4.;
  for (auto p : S3Dinverted) {if (p.z() > maxZ) maxZ = p.z();}
  std::cout << "maxZ: " << maxZ << std::endl;

  // Initialize bbox and number of cells blocks corresponding to resolution
  const std::vector<std::vector<double>> bbox = make_bbox(pointcloud);
  const std::vector<int> CELLROWSCOLS = make_cells(bbox, resolution);
  const int CELLROWS = CELLROWSCOLS[0];
  const int CELLCOLS = CELLROWSCOLS[1];
  std::cout << "CELLROWS: " << CELLROWS << "\n";
  std::cout << "CELLCOLS: " << CELLCOLS << "\n";

  const double offsetX = bbox[0][0];
  const double offsetY = bbox[0][1];

  // Initialize cloth as regularly distributed grid of particles
  std::vector<std::vector<Point>> cloth(CELLROWS,std::vector<Point>(CELLCOLS));

  std::cout << "cloth size: " << cloth.size() << std::endl;

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

  // Calculate particle elevation parameters
  std::vector<std::vector<std::vector<double>>> zParams(CELLROWS,std::vector<std::vector<double>>(CELLCOLS));
  //std::vector<std::vector<double>> zParams; // Initialise vector to store z params

  // Construct and query kd-tree:
  // https://doc.cgal.org/latest/Spatial_searching/index.html#title5
  Tree tree(S2D.begin(), S2D.end());
  const unsigned int N = 1;
  std::vector<Point>::iterator itS2D;
  //std::vector<double> displacement = {0., 0., 4.};
  double disp = 0.5;
  int actual_cloth_size = 0;
  for (int i=0; i < CELLROWS; i++) {
      for (int j=0; j < CELLCOLS; j++) {
          actual_cloth_size++;
          double pZmin, pZprev, pZcur;
          std::vector<double> zParam;
          //zParam[0] = pZmin; zParam[1] = pZprev; zParam[2] = pZcur;

          //std::cout << "i, j:               " << i << ", " << j << std::endl;
          //std::cout << "cloth[i][j]:        " << cloth[i][j] << std::endl;
          Point query_point = Point(cloth[i][j].x(), cloth[i][j].y(), 0.);
          Neighbor_search search_result(tree, query_point, N);
          for (auto res: search_result) {
              Point neighbour_point = res.first;
              //double distance = res.second;
              //double distanceSquare = std::sqrt(res.second);
              //std::cout << "query_point:        " << query_point << std::endl;
              //std::cout << "neighbour_point:    " << neighbour_point << std::endl;
              //std::cout << "neighbour_point z:    " << neighbour_point.z() << std::endl;

              itS2D = std::find(S2D.begin(), S2D.end(), neighbour_point);
              int idxS3D = std::distance(S2D.begin(), itS2D);
              pZmin = S3Dinverted[idxS3D].z();
              zParams[i][j].push_back(pZmin);
              pZcur = z0;
              pZprev = z0 + disp;
              zParams[i][j].push_back(pZprev);
              zParams[i][j].push_back(pZcur);
              zParams[i][j].push_back(1.); // Insert movable variable
              //std::cout << "Sinverse p:         " << S3Dinverted[idxS3D] << std::endl;
              //pZmin = neighbour_point.z(); zParam.push_back(pZmin);
              //std::cout << "distance:           " << distance << std::endl;
              //std::cout << "distanceSquare:     " << distanceSquare << std::endl;
          }
      }
  }
    int cone=0; int ctwo=0; int cthree=0; int cfour=0; int cfive=0; int csix=0; int cseven=0; int ceight=0; int cnine=0;
  //int movable = 0;
  double deltaZ = disp;
  int deltaZstopper = 0;
  while (deltaZ > epsilon_zmax) {

      // external forces, apply to all movable p
      for (int i = 0; i < CELLROWS; i++) {
          for (int j = 0; j < CELLCOLS; j++) {
              if (zParams[i][j][3] == 1.) {
                  //std::cout << "[i][j]:         [" << i << "][" << j << "]" << std::endl;
                  //std::cout << "zParams[i][j]:    (" << zParams[i][j][0] << ", " << zParams[i][j][1] << ", "
                  //          << zParams[i][j][2] << ")" << std::endl;
                  //std::cout << "after one iteration, the result of step 10-13 is:" << std::endl;
                  double pzmin = zParams[i][j][0];
                  double pzprev = zParams[i][j][1];
                  double pzcur = zParams[i][j][2];

                  double tmp = pzcur;
                  pzcur = (pzcur - pzprev) + pzcur;
                  pzprev = tmp;
                  if (pzcur <= pzmin) {pzcur = pzmin; zParams[i][j][3] = 0.;}
                  zParams[i][j][1] = pzprev;
                  zParams[i][j][2] = pzcur;
                  //std::cout << "pzcur:      " << pzcur << std::endl;
                  //std::cout << "pzprev:     " << pzprev << std::endl;
                  //if (zParams[i][j][3] == 0.) { std::cout << "Indeed, this one is movable...\n"; movable++;}
              }
          }
      }

      // internal forces, process once each set neighbours e of adjacent particles

      for (int i = 0; i < CELLROWS; i++) {
          for (int j = 0; j < CELLCOLS; j++) {
              //std::cout << "[i][j]:         [" << i << "][" << j << "]" << std::endl;
              if (zParams[i][j][3] == 1.) {
                  if ((i > 0) && (i < (CELLROWS - 1)) && (j > 0) && (j < (CELLCOLS - 1))) {
                      cfive++;
                      std::vector<std::vector<double>> ni;
                      std::vector<double> n1, n2, n3, n4;
                      n1 = zParams[i-1][j]; ni.push_back(n1);
                      n2 = zParams[i+1][j]; ni.push_back(n2);
                      n3 = zParams[i][j-1]; ni.push_back(n3);
                      n4 = zParams[i][j+1]; ni.push_back(n4);
                      // n1
                      if (zParams[i-1][j][3] == 0.){
                          if (zParams[i-1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                          }
                          else if (zParams[i-1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                          }
                      }
                      if (zParams[i-1][j][3] == 1.){
                          if (zParams[i-1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                              zParams[i-1][j][2] += disp;
                          }
                          else if (zParams[i-1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                              zParams[i-1][j][2] -= disp;
                          }
                      }
                      // n2
                      if (zParams[i+1][j][3] == 0.){
                          if (zParams[i+1][j][2] < cloth[i][j].z()){
                          cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                        }
                          else if (zParams[i+1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                          }
                      }
                      if (zParams[i+1][j][3] == 1.){
                          if (zParams[i+1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                              zParams[i+1][j][2] += disp;
                          }
                          else if (zParams[i+1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                              zParams[i+1][j][2] -= disp;
                          }
                      }
                      // n3
                      if (zParams[i][j-1][3] == 0.){
                          if (zParams[i][j-1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                          }
                          else if (zParams[i][j-1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                          }
                      }
                      if (zParams[i][j-1][3] == 1.){
                          if (zParams[i][j-1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                              zParams[i][j-1][2] += disp;
                          }
                          else if (zParams[i][j-1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                              zParams[i][j-1][2] -= disp;
                          }
                      }
                      //n4
                      if (zParams[i][j+1][3] == 0.){
                          if (zParams[i][j+1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                          }
                          else if (zParams[i][j+1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                          }
                      }
                      if (zParams[i][j+1][3] == 1.){
                          if (zParams[i][j+1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                              zParams[i][j+1][2] += disp;
                          }
                          else if (zParams[i][j+1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                              zParams[i][j+1][2] -= disp;
                          }
                      }

                  } // 5
                  else if ((i == 0) && (j > 0) && (j < (CELLCOLS - 1))) {
                      ctwo++;
                      std::vector<std::vector<double>> ni;
                      std::vector<double> n1, n2, n3;
                      n1 = zParams[i][j-1]; ni.push_back(n1);
                      n2 = zParams[i][j+1]; ni.push_back(n2);
                      n3 = zParams[i+1][j]; ni.push_back(n3);
                      // n1
                      if (zParams[i][j-1][3] == 0.){
                          if (zParams[i][j-1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                          }
                          else if (zParams[i][j-1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                          }
                      }
                      if (zParams[i][j-1][3] == 1.){
                          if (zParams[i][j-1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                              zParams[i][j-1][2] += disp;
                          }
                          else if (zParams[i][j-1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                              zParams[i][j-1][2] -= disp;
                          }
                      }
                      // n2
                      if (zParams[i][j+1][3] == 0.){
                          if (zParams[i][j+1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                          }
                          else if (zParams[i][j+1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                          }
                      }
                      if (zParams[i][j+1][3] == 1.){
                          if (zParams[i][j+1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                              zParams[i][j+1][2] += disp;
                          }
                          else if (zParams[i][j+1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                              zParams[i][j+1][2] -= disp;
                          }
                      }
                      // n3
                      if (zParams[i+1][j][3] == 0.){
                          if (zParams[i+1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                          }
                          else if (zParams[i+1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                          }
                      }
                      if (zParams[i+1][j][3] == 1.){
                          if (zParams[i+1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                              zParams[i+1][j][2] += disp;
                          }
                          else if (zParams[i+1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                              zParams[i+1][j][2] -= disp;
                          }
                      }

                  } // 2
                  else if ((i > 0) && (i < (CELLROWS - 1)) && (j == 0)) {
                      cfour++;
                      std::vector<std::vector<double>> ni;
                      std::vector<double> n1, n2, n3;
                      n1 = zParams[i][j+1]; ni.push_back(n1);
                      n2 = zParams[i-1][j]; ni.push_back(n2);
                      n3 = zParams[i+1][j]; ni.push_back(n3);
                      // n1
                      if (zParams[i][j+1][3] == 0.){
                          if (zParams[i][j+1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                          }
                          else if (zParams[i][j+1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                          }
                      }
                      if (zParams[i][j+1][3] == 1.){
                          if (zParams[i][j+1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                              zParams[i][j+1][2] += disp;
                          }
                          else if (zParams[i][j+1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                              zParams[i][j+1][2] -= disp;
                          }
                      }
                      // n2
                      if (zParams[i-1][j][3] == 0.){
                          if (zParams[i-1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                          }
                          else if (zParams[i-1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                          }
                      }
                      if (zParams[i-1][j][3] == 1.){
                          if (zParams[i-1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                              zParams[i-1][j][2] += disp;
                          }
                          else if (zParams[i-1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                              zParams[i-1][j][2] -= disp;
                          }
                      }
                      // n3
                      if (zParams[i+1][j][3] == 0.){
                          if (zParams[i+1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                          }
                          else if (zParams[i+1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                          }
                      }
                      if (zParams[i+1][j][3] == 1.){
                          if (zParams[i+1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                              zParams[i+1][j][2] += disp;
                          }
                          else if (zParams[i+1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                              zParams[i+1][j][2] -= disp;
                          }
                      }

                  } // 4
                  else if ((i > 0) && (i < (CELLROWS - 1)) && (j == (CELLCOLS - 1))) {
                      csix++;
                      std::vector<std::vector<double>> ni;
                      std::vector<double> n1, n2, n3;
                      n1 = zParams[i][j-1]; ni.push_back(n1);
                      n2 = zParams[i-1][j]; ni.push_back(n2);
                      n3 = zParams[i+1][j]; ni.push_back(n3);
                      // n1
                      if (zParams[i][j-1][3] == 0.){
                          if (zParams[i][j-1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                          }
                          else if (zParams[i][j-1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                          }
                      }
                      if (zParams[i][j-1][3] == 1.){
                          if (zParams[i][j-1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                              zParams[i][j-1][2] += disp;
                          }
                          else if (zParams[i][j-1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                              zParams[i][j-1][2] -= disp;
                          }
                      }
                      // n2
                      if (zParams[i-1][j][3] == 0.){
                          if (zParams[i-1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                          }
                          else if (zParams[i-1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                          }
                      }
                      if (zParams[i-1][j][3] == 1.){
                          if (zParams[i-1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                              zParams[i-1][j][2] += disp;
                          }
                          else if (zParams[i-1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                              zParams[i-1][j][2] -= disp;
                          }
                      }
                      //n3
                      if (zParams[i+1][j][3] == 0.){
                          if (zParams[i+1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                          }
                          else if (zParams[i+1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                          }
                      }
                      if (zParams[i+1][j][3] == 1.){
                          if (zParams[i+1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                              zParams[i+1][j][2] += disp;
                          }
                          else if (zParams[i+1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                              zParams[i+1][j][2] -= disp;
                          }
                      }

                  } // 6
                  else if ((i == (CELLROWS - 1)) && (j > 0) && (j < (CELLCOLS - 1))) {
                      ceight++;
                      std::vector<std::vector<double>> ni;
                      std::vector<double> n1, n2, n3;
                      n1 = zParams[i][j-1]; ni.push_back(n1);
                      n2 = zParams[i][j+1]; ni.push_back(n2);
                      n3 = zParams[i-1][j]; ni.push_back(n3);
                      // n1
                      if (zParams[i][j-1][3] == 0.){
                          if (zParams[i][j-1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                          }
                          else if (zParams[i][j-1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                          }
                      }
                      if (zParams[i][j-1][3] == 1.){
                          if (zParams[i][j-1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                              zParams[i][j-1][2] += disp;
                          }
                          else if (zParams[i][j-1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                              zParams[i][j-1][2] -= disp;
                          }
                      }
                      // n2
                      if (zParams[i][j+1][3] == 0.){
                          if (zParams[i][j+1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                          }
                          else if (zParams[i][j+1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                          }
                      }
                      if (zParams[i][j+1][3] == 1.){
                          if (zParams[i][j+1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                              zParams[i][j+1][2] += disp;
                          }
                          else if (zParams[i][j+1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                              zParams[i][j+1][2] -= disp;
                          }
                      }
                      // n3
                      if (zParams[i-1][j][3] == 0.){
                          if (zParams[i-1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                          }
                          else if (zParams[i-1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                          }
                      }
                      if (zParams[i-1][j][3] == 1.){
                          if (zParams[i-1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                              zParams[i-1][j][2] += disp;
                          }
                          else if (zParams[i-1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                              zParams[i-1][j][2] -= disp;
                          }
                      }
                  } // 8
                  else if ((i == 0) && (j == 0)) {
                      cone++;
                      std::vector<std::vector<double>> ni;
                      std::vector<double> n1, n2;
                      n1 = zParams[i][j+1]; ni.push_back(n1);
                      n2 = zParams[i+1][j]; ni.push_back(n2);
                      // n1
                      if (zParams[i][j+1][3] == 0.){
                          if (zParams[i][j+1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                          }
                          else if (zParams[i][j+1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                          }
                      }
                      if (zParams[i][j+1][3] == 1.){
                          if (zParams[i][j+1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                              zParams[i][j+1][2] += disp;
                          }
                          else if (zParams[i][j+1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                              zParams[i][j+1][2] -= disp;
                          }
                      }
                      // n2
                      if (zParams[i+1][j][3] == 0.){
                          if (zParams[i+1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                          }
                          else if (zParams[i+1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                          }
                      }
                      if (zParams[i+1][j][3] == 1.){
                          if (zParams[i+1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                              zParams[i+1][j][2] += disp;
                          }
                          else if (zParams[i+1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                              zParams[i+1][j][2] -= disp;
                          }
                      }

                  } // 1
                  else if ((i == 0) && (j == (CELLCOLS - 1))) {
                      cthree++;
                      std::vector<std::vector<double>> ni;
                      std::vector<double> n1, n2;
                      n1 = zParams[i][j-1]; ni.push_back(n1);
                      n2 = zParams[i+1][j]; ni.push_back(n2);
                      // n1
                      if (zParams[i][j-1][3] == 0.){
                          if (zParams[i][j-1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                          }
                          else if (zParams[i][j-1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                          }
                      }
                      if (zParams[i][j-1][3] == 1.){
                          if (zParams[i][j-1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                              zParams[i][j-1][2] += disp;
                          }
                          else if (zParams[i][j-1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                              zParams[i][j-1][2] -= disp;
                          }
                      }
                      // n2
                      if (zParams[i+1][j][3] == 0.){
                          if (zParams[i+1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                          }
                          else if (zParams[i+1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                          }
                      }
                      if (zParams[i+1][j][3] == 1.){
                          if (zParams[i+1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                              zParams[i+1][j][2] += disp;
                          }
                          else if (zParams[i+1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                              zParams[i+1][j][2] -= disp;
                          }
                      }

                  } // 3
                  else if ((i == (CELLROWS - 1)) && (j == 0)) {
                      cseven++;
                      std::vector<std::vector<double>> ni;
                      std::vector<double> n1, n2;
                      n1 = zParams[i-1][j]; ni.push_back(n1);
                      n2 = zParams[i][j+1]; ni.push_back(n2);
                      // n1
                      if (zParams[i-1][j][3] == 0.){
                          if (zParams[i-1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                          }
                          else if (zParams[i-1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                          }
                      }
                      if (zParams[i-1][j][3] == 1.){
                          if (zParams[i-1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                              zParams[i-1][j][2] += disp;
                          }
                          else if (zParams[i-1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                              zParams[i-1][j][2] -= disp;
                          }
                      }
                      // n2
                      if (zParams[i][j+1][3] == 0.){
                          if (zParams[i][j+1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                          }
                          else if (zParams[i][j+1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                          }
                      }
                      if (zParams[i][j+1][3] == 1.){
                          if (zParams[i][j+1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                              zParams[i][j+1][2] += disp;
                          }
                          else if (zParams[i][j+1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                              zParams[i][j+1][2] -= disp;
                          }
                      }

                  } // 7
                  else if ((i == (CELLROWS - 1)) && (j == (CELLCOLS - 1))) {
                      cnine++;
                      std::vector<std::vector<double>> ni;
                      std::vector<double> n1, n2;
                      n1 = zParams[i][j-1]; ni.push_back(n1);
                      n2 = zParams[i-1][j]; ni.push_back(n2);
                      // n1
                      if (zParams[i][j-1][3] == 0.){
                          if (zParams[i][j-1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                          }
                          else if (zParams[i][j-1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                          }
                      }
                      if (zParams[i][j-1][3] == 1.){
                          if (zParams[i][j-1][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                              zParams[i][j-1][2] += disp;
                          }
                          else if (zParams[i][j-1][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                              zParams[i][j-1][2] -= disp;
                          }
                      }
                      // n2
                      if (zParams[i-1][j][3] == 0.){
                          if (zParams[i-1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                          }
                          else if (zParams[i-1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                          }
                      }
                      if (zParams[i-1][j][3] == 1.){
                          if (zParams[i-1][j][2] < cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() - disp);
                              zParams[i-1][j][2] += disp;
                          }
                          else if (zParams[i-1][j][2] > cloth[i][j].z()){
                              cloth[i][j] = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z() + disp);
                              zParams[i-1][j][2] -= disp;
                          }
                      }

                  } // 9
              }
          }
      }
      // Calculate the max deltaZ
      for (int i=0; i < CELLROWS; i++) {
          for (int j=0; j < CELLCOLS; j++) {
              //std::cout << "(zParams[i][j][2] - zParams[i][j][1]):      " << "(" << zParams[i][j][2] << " " << "- " << zParams[i][j][1] << ")\n";
              if ((zParams[i][j][2] - zParams[i][j][1]) > deltaZ) {

                  //std::cout << "(zParams[i][j][2] - zParams[i][j][1]):      " << "(" << zParams[i][j][2] << " " << "- " << zParams[i][j][1] << ")\n";
                  deltaZ = zParams[i][j][2] - zParams[i][j][1];
              }
              if ((zParams[i][j][2] - zParams[i][j][1]) == deltaZ){deltaZstopper++;}
          }
      }

      if (deltaZstopper >= 500){
          std::cout << "STOP!" << "\n";
          break;}

  }

  std::vector<Point> clothVec;
    for (int i=0; i < CELLROWS; i++) {
        for (int j=0; j < CELLCOLS; j++) {
            clothVec.push_back(cloth[i][j]);
            // if (i == 0) {std::cout << "cloth[0][j]:     " << cloth[i][j] << "\n";}
            // if (j == 0) {std::cout << "cloth[i][0]:     " << cloth[i][j] << "\n";}

        }
    }

    for (int i=0; i < CELLROWS; i++) {
        for (int j=0; j < CELLCOLS; j++) {

            // if (i == (CELLROWS - 1)) {std::cout << "cloth[max][j]:      " << cloth[i][j] << "\n";}
            // if (j == (CELLCOLS - 1)) {std::cout << "cloth[i][max]:      " << cloth[i][j] << "\n";}
        }
    }

  std::vector<int> class_labels;
  Tree tree3d(clothVec.begin(), clothVec.end());
  int idx=0;
  for (auto p : S3Dinverted) {
      Point query_point = p;
      Neighbor_search search_result(tree3d, query_point, N);
      for (auto res : search_result) {
          Point neighbour_point = res.first;
          double distance = res.second;
          if (distance < (epsilon_ground*epsilon_ground)) {class_labels.push_back(2);
              //std::cout << "groundpoint found...\n";
          }
          else {class_labels.push_back(1);}
      }
      idx++;
  }

/*
    std::vector<Point>::iterator itS3D;
    for (int i=0; i < CELLROWS; i++) {
        for (int j=0; j < CELLCOLS; j++) {
            Point query_point = Point(cloth[i][j].x(), cloth[i][j].y(), cloth[i][j].z());
            Neighbor_search search_result(tree3d, query_point, N);
            for (auto res: search_result) {
                Point neighbour_point = res.first;
                double distance = std::sqrt(res.second);
                if (distance < epsilon_ground) {

                    itS3D = std::find(S3Dinverted.begin(), S3Dinverted.end(), neighbour_point);
                    int idxS3D = std::distance(S3Dinverted.begin(), itS3D);

                }
            }
        }
    }
    */
  std::cout << "all the counts: \n one: " << cone << "\n two: " << ctwo << "\n three: " << cthree << "\n four: " << cfour << "\n five: ";
  std::cout << cfive << "\n six: " << csix << "\n seven: " << cseven << "\n eight: " << ceight << "\n nine: " << cnine << "\n";




  std::cout << "actual cloth size: " << actual_cloth_size << std::endl;
  //std::cout << "movables: " << movable << std::endl;

  double pZmin, pZprev, pZcur;
  std::vector<double> zParam;
  pZmin = 0., pZprev = 1., pZcur = 2.;
  zParam.push_back(pZmin); zParam.push_back(pZprev); zParam.push_back(pZcur);

  std::cout << "zParam vector:      " << zParam[0] << " " << zParam[1] << " " << zParam[2] << " " << std::endl;
  std::vector<double> difference;
  //std::transform(zParam.begin(), zParam.end(),
  //                             displacement.begin(), zParam.begin(),
  //                             std::minus<double>());

    std::cout << "zParam contains:";
    for (std::vector<double>::iterator it=zParam.begin(); it!=zParam.end(); ++it)
        std::cout << ' ' << *it;
    std::cout << '\n';

  //  displace_pt(cloth[0][0], displacement);

    std::cout << "cloth[0][0] contains:";
    for (std::vector<double>::iterator it=zParam.begin(); it!=zParam.end(); ++it)
        std::cout << ' ' << *it;
    std::cout << '\n';

  //-- TIP
  //-- write the results to a new LAS file

  write_lasfile(jparams["output_las"], pointcloud, class_labels);
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

/*
 * Notes and such
 *
 *
 *     int count = 0;

    for (auto label : class_labels) {
        //std::cout << label << std::endl;
        if (label == 2) {
            count++;
            //std::cout << label << std::endl;
        }
    }
    std::cout << "the occurrence of label 2:        " << count << std::endl;
 *
  DT::Face_handle triangle = dt.locate(Point(5,5,0));
  //-- get the 3 vertices of the triangle:
  //DT::Vertex_handle v0 = triangle->vertex(0);
  //DT::Vertex_handle v1 = triangle->vertex(1);
  //DT::Vertex_handle v2 = triangle->vertex(2);
  // get the coordinates of the three vertices:
  // std::cout << "v0 has the coordinates ( " << v0->point().x() << "  " << v0->point().y() << " " << v0->point().z() << " )" << std::endl;
  // std::cout << "v1 has the coordinates ( " << v1->point().x() << "  " << v1->point().y() << " " << v1->point().z() << " )" << std::endl;
  // std::cout << "v2 has the coordinates ( " << v2->point().x() << "  " << v2->point().y() << " " << v2->point().z() << " )" << std::endl;

  //-- TIP CGAL compute squared distance between two points: [https://doc.cgal.org/latest/Kernel_23/group__squared__distance__grp.html#ga1ff73525660a052564d33fbdd61a4f71]
  //std::cout << "the squared distance between v0 and v1 is: " << CGAL::squared_distance(v0->point(), v1->point()) << std::endl;

 *
 *     int count = 0;
    for (size_t row=0; row < rows; row += resolution) {
        for (size_t col=0; col < cols; col += resolution) {
            //std::cout << "zeh each row and cols: " << row << ", " << col << std::endl;
            count++;
            std::vector<double> localMin;
            for (size_t i=0; i < resolution; ++i) {
                for (size_t j=0; j < resolution; ++j) {
                    //std::cout << "point: " << row + i << ", " << col + j << std::endl;
                    localMin.push_back(pointcloud[row + i][col + j]);
                }
            }
            std::vector<double>::iterator result = std::min_element(localMin.begin(), localMin.end());
            //std::cout << "localMin: " << result[0] << std::endl;
        }
    }
    std::cout << "the count: " << count << std::endl;
 *
 *   DT dtT;
  dtT.insert(Point(0,0,4));
  dtT.insert(Point(10,0,3));
  dtT.insert(Point(0,10,0));
  //-- Find triangle that intersects a given point: [https://doc.cgal.org/latest/Triangulation_2/classCGAL_1_1Triangulation__2.html#a940567120751e7864c7b345eaf756642]

  Point p1, p2, p0, p4, p5, p6;
  p0 = Point(0,0,0);
  p1 = Point(4, 8, 10);
  p2 = Point(9, 2, 7);
  p4 = Point(0.,0.,0.);
  p5 = Point(1.,1.,1.);
  p6 = Point(2.,2.,2.);

  Vector a, b;
  a = Vector(p0, p1);
  b = Vector(p0, p2);

  std::vector<Point> test;
  test.push_back(p4);
  test.push_back(p5);
  test.push_back(p6);

  std::vector<Point>::iterator itt;
  itt = std::find(test.begin(), test.end(), p5);
    if (itt != test.end())
        std::cout << "Element Found" << std::endl;
    else
        std::cout << "Element Not Found" << std::endl;

    int index = std::distance(test.begin(), itt);
    std::cout << "and test index is: " << index << std::endl;

 *
 *   std::vector<Point>::const_iterator it;
  it = std::find(pointcloud.begin(), pointcloud.end(), p4);

    if (it != test.end())
        std::cout << "Element Found" << std::endl;
    else
        std::cout << "Element Not Found" << std::endl;

    int index = std::distance(test.begin(), it);
    std::cout << "and its index is: " << index << std::endl;

    int indexT = 4;
    Point test_idx;
    test_idx = Point (0.,0.,0.);
    if (test[4] == Point (0.,0.,0.)) {std::cout << "testing invalid index WORKS " << std::endl;}

 *   int count = 0;
  int countWeird = 0;
  int countOne = 0;

  for (auto label : class_labels) {
      //std::cout << label << std::endl;
      if (label == 2) {
          count++;
          //std::cout << label << std::endl;
      }
      if (label != 1) {countWeird++; std::cout << label;}
      if (label == 1) {countOne++;}
  }
  std::cout << "the occurrence of label 2:        " << count << std::endl;
  std::cout << "the occurrence of non-label 1:        " << countWeird << std::endl;
  std::cout << "the occurrence of label 1:        " << countOne << std::endl;
  std::cout << "the actual occurrence of label 2: " << 19 * 15 << std::endl;
  int dot;

  dot = a * b;
  std::cout << "the inner/dot product of a * b = " << dot << std::endl;

  Point c1, c2, c0, c4;
  c0 = Point(0,0,0);
  c1 = Point(2, 3, 4);
  c2 = Point(5, 6, 7);
  Vector ca, cb, cross, sub;
  ca = Vector(c0, c1);
  cb = Vector(c0, c2);
  sub = c2 - c1;
  c4 = Point(sub[0], sub[1], sub[2]);
  double dot2;
  dot2 = sub * ca;

  std::cout << "subtraction vector of points c2 - c1 = " << sub << " which is the same as: " << c4 << std::endl;
  std::cout << "and the dot product of the result * ca = (" << sub << ") * " << ca << " = " << dot2 << std::endl;

  cross = CGAL::cross_product(ca, cb);
  std::cout << "the cross product of a x b = " << cross << std::endl;



 */