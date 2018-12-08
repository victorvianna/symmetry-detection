#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/adjacency_list.h>
#include <vector>
#include "mean_shift.h"
#include "kernel.h"

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

int main(int argc, char *argv[])
{
  //test_mesh_display();
  test_clustering();
  return 0;
}
