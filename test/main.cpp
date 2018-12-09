#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/adjacency_list.h>
#include <vector>
#include "mean_shift.h"
#include "kernel.h"
#include "flann/flann.hpp"

void test_mesh_display(){
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;

  std::vector<std::vector<int>> A;

  // Load a mesh in OFF format
  igl::readOFF("mesh/bunny.off", V, F);

  igl::adjacency_list(F, A);

  int edges = 0;
  for (auto vertex : A){
    edges += vertex.size();
  }
  std::cout << "Edges from adjacency = " << edges / 2 << std::endl;

  std::cout << "Edges from faces = " << (F.rows() * 3) / 2 << std::endl;

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.launch();

  igl::adjacency_list(F, A);

}

void test_clustering(){
    using namespace std;
    MeanShift *msp = new MeanShift(epanechnikov_kernel);
    double kernel_bandwidth = 3;

    vector<vector<double> > points = {{0, 0, 0}, {1, 1, 1}};
    vector<Cluster> clusters = msp->cluster(points, kernel_bandwidth);
}

void test_nearest_neighbors(){
    using namespace flann;

    double* dataset_ = new double[8]{
            0, 0,
            1, 1,
            2, 2,
            3, 3
    };
    Matrix<double> dataset(dataset_, 4, 2);
    double* query_ = new double[2]{
            0.5, 0.5
    };
    Matrix<double> query(query_, 1, 2);

    std::vector< std::vector<int>> indices;
    std::vector<std::vector<double>> dists;

    // construct an randomized kd-tree index using 4 kd-trees

    Index<L2<double>> index(dataset, flann::KDTreeIndexParams(4));
    index.buildIndex();

    // search around a certain radius
    double radius = 0.6;
    index.radiusSearch(query, indices, dists, radius, flann::SearchParams(128));

    delete[] dataset.ptr();
    delete[] query.ptr();
}

int main(int argc, char *argv[])
{
  //test_mesh_display();
  test_clustering();
  test_nearest_neighbors();
  return 0;
}
