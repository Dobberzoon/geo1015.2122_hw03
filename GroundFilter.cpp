/*
  GEO1015.2021
  hw03 
  --
  [YOUR NAME] 
  [YOUR STUDENT NUMBER] 
  [YOUR NAME] 
  [YOUR STUDENT NUMBER] 
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

// -- CGAL other
typedef CGAL::Search_traits_3<Kernel> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
typedef Neighbor_search::Tree Tree;


// TODO
/*
 * This seems overcomplicated, try to keep things as simple as possible. The simplest (and fastest) approach would be to
1) define the virtual grid where each cell can hold one point.
2) iterate over the pointcloud vector and determine for each point using it's x and y coordinates in which gridcell it falls (this is a simple one line formula). Then either put the point in that gridcell if it is empty or overwrite the point in the gridcell if the current point is lower (compare using the z coordinate).
3) Finally with the points that are now in the virtual grid build the initial TIN.

This should be very fast (time complexity O(N)) and no KD tree is needed.
 */

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

  for (auto p: pointcloud) {if (p.x() == 0. || p.y() == 0. || p.z() == 0.) {std::cout << "outlier detected: " << p << std::endl;}}

  double resolution = jparams["resolution"];
  double distance = jparams["distance"];
  double angle = jparams["angle"];
  std::string output_las = jparams["output_las"];

  std::cout << "output_las: " << output_las << std::endl;
  std::cout << "resolution: " << resolution << std::endl;
  std::cout << "distance: " << distance << std::endl;
  std::cout << "angle: " << angle << std::endl;

  // Initialize bbox and number of cells blocks corresponding to resolution
  const std::vector<std::vector<double>> bbox = make_bbox(pointcloud);
  //std::cout << "the bbox: " << "upperLeft: (" << bbox[0][0] << ", " << bbox[0][1] << ")" << " lowerRight: (" << bbox[1][0] << ", " << bbox[1][1] << ")" << std::endl;
  const std::vector<int> CELLROWSCOLS = make_cells(bbox, resolution);
  std::cout << "the number of CELLROWS: " << CELLROWSCOLS[0] << ", CELLCOLS: " << CELLROWSCOLS[1] << std::endl;
  const int CELLROWS = CELLROWSCOLS[0];
  const int CELLCOLS = CELLROWSCOLS[1];
  //const int rows = std::ceil(bbox[1][1] - bbox[0][1]);
  //const int cols = std::ceil(bbox[1][0] - bbox[0][0]);

  const double offsetX = bbox[0][0];
  const double offsetY = bbox[0][1];

  //std::cout << "rows, cols: " << rows << ", " << cols << std::endl;
  std::cout << "offset x, y: " << offsetX << ", " << offsetY << std::endl;

  // Initialize virtual grid to store initial ground points
  std::vector<std::vector<Point>> vGrid(CELLROWS,std::vector<Point>(CELLCOLS));

  // The virtual grid will have one point per cellblock, which will contain lowest elevation value.
  for (auto p : pointcloud) {
      int cellX, cellY;
      cellX = std::round((p.y() - offsetY) / resolution);
      cellY = std::round((p.x() - offsetX) / resolution);
      //std::cout << "x, y: " << cellX << ", " << cellY << std::endl;
      //std::cout << "what's inside this cell: " << vGrid[cellX][cellY] << std::endl;

      // This pointcloud dataset contains noise/outliers that have some or all p(x, y, z) = ~(0., 0., 0.)
      // the optional filter  && (p[0] != 0. && p[1] != 0. && p[2] != 0.) or  && (p[0] != 0. || p[1] != 0. || p[2] != 0.)
      if ((vGrid[cellX][cellY][0] == 0. && vGrid[cellX][cellY][1] == 0. && vGrid[cellX][cellY][2] == 0.)) {vGrid[cellX][cellY] = p;}
      else {
          if (p[2] < vGrid[cellX][cellY][2]) {vGrid[cellX][cellY] = p;}
      }

  }

  //-- TIP CGAL triangulation -> https://doc.cgal.org/latest/Triangulation_2/index.html
  //-- Insert points in a triangulation: [https://doc.cgal.org/latest/Triangulation_2/classCGAL_1_1Triangulation__2.html#a1025cd7e7226ccb44d82f0fb1d63ad4e]
  DT dt;

  // Insertion of initial points from vGrid into DT
  for (int i=0; i < CELLROWS; i++) {
      for (int j=0; j < CELLCOLS; j++) {
          //std::cout << "vGrid[" << i << "][" << j << "] = " << vGrid[i][j] << std::endl;
          dt.insert(vGrid[i][j]);
      }
  }



    /*
    int count = 0;
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
    */

  DT dtT;
  dtT.insert(Point(0,0,4));
  dtT.insert(Point(10,0,3));
  dtT.insert(Point(0,10,0));
  //-- Find triangle that intersects a given point: [https://doc.cgal.org/latest/Triangulation_2/classCGAL_1_1Triangulation__2.html#a940567120751e7864c7b345eaf756642]
  DT::Face_handle triangle = dtT.locate(Point(5,5,1000));
  //-- get the 3 vertices of the triangle: 
  DT::Vertex_handle v0 = triangle->vertex(0);
  DT::Vertex_handle v1 = triangle->vertex(1);
  DT::Vertex_handle v2 = triangle->vertex(2);
  // get the coordinates of the three vertices:
  std::cout << "v0 has the coordinates ( " << v0->point().x() << "  " << v0->point().y() << " " << v0->point().z() << " )" << std::endl;
  std::cout << "v1 has the coordinates ( " << v1->point().x() << "  " << v1->point().y() << " " << v1->point().z() << " )" << std::endl;
  std::cout << "v2 has the coordinates ( " << v2->point().x() << "  " << v2->point().y() << " " << v2->point().z() << " )" << std::endl;

  //-- TIP CGAL compute squared distance between two points: [https://doc.cgal.org/latest/Kernel_23/group__squared__distance__grp.html#ga1ff73525660a052564d33fbdd61a4f71]
  std::cout << "the squared distance between v0 and v1 is: " << CGAL::squared_distance(v0->point(), v1->point()) << std::endl;
  
  //-- TIP
  //-- write the results to a new LAS file
  //std::vector<int> class_labels;

  //write_lasfile(jparams["output_las"], pointcloud, class_labels);
}



/*
  bbox part_wassenaar_ahn3.laz
  min x y z:                  84508.948 460678.587 4.881
  max x y z:                  84599.998 460749.999 33.168

  expected number of cells with resolution=5
  CELLCOLS -> 84600 - 84509 = 91 -> 91 / 5 = 18.2 -> 19
  CELLSROWS -> 460750 - 460679 = 71 -> 71 / 5 = 14.2 -> 15
 */




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


  // double resolution = jparams["resolution"];
  // double epsilon_zmax = jparams["epsilon_zmax"];
  // double epsilon_ground = jparams["epsilon_ground"];
  // std::string output_las = jparams["output_las"];

  // //-- print the first 5 points in the pointcloud, which are CGAL Point_3
  // //-- https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Point__3.html
  // int i = 0;
  // for (auto p : pointcloud) {
  //   std::cout << "(" << p.x() << ", " << p.y() << ", " << p.z()  << ")" << std::endl;
  //   i++;
  //   if (i == 5)
  //     break;
  // }

  //-- TIP
  //-- construct and query kd-tree:
  // https://doc.cgal.org/latest/Spatial_searching/index.html#title5
  //
  // Tree tree(pointcloud.begin(), pointcloud.end());
  // const unsigned int N = 1;
  // Point query_point = Point(0,0,0);
  // Neighbor_search search_result(tree, query_point, N);
  //
  // for(auto res : search_result) {
  //   Point neighbour_point = res.first;
  //   double distance = res.second;
  // }

  //-- TIP
  //-- write the results to a new LAS file
  // std::vector<int> class_labels;
  // write_lasfile(jparams["output_las"], pointcloud, class_labels);
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

