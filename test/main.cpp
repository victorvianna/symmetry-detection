#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/adjacency_list.h>
#include <vector>

Eigen::MatrixXd V;
Eigen::MatrixXi F;

std::vector<std::vector<int>> A;

int main(int argc, char *argv[])
{
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

  return 0;
}
