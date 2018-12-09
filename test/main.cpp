#include <igl/avg_edge_length.h>
#include <igl/parula.h>
#include <igl/principal_curvature.h>
#include <igl/read_triangle_mesh.h>
#include <igl/adjacency_list.h>
#include <igl/opengl/glfw/Viewer.h>

Eigen::MatrixXd vertices;
Eigen::MatrixXi faces;

Eigen::MatrixXd minCurvDir, maxCurvDir;
Eigen::VectorXd minCurvVal, maxCurvVal;

//std::vector<std::vector<int>> adjacencyList;

const Eigen::RowVector3d red(0.8, 0.2, 0.2), blue(0.2, 0.2, 0.8);
double length;

// Remove the segment lines added to the mesh
void clear_lines(igl::opengl::glfw::Viewer& viewer){
    viewer.data().lines.resize(0, 9);
}

// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier){
    switch (key){
        case '1':
            clear_lines(viewer);
            // Draw a blue segment parallel to the minimal curvature direction
            viewer.data().add_edges(vertices + minCurvDir * length, vertices - minCurvDir * length, blue);
            return true;
        case '2':
            clear_lines(viewer);
            // Draw a red segment parallel to the maximal curvature direction
            viewer.data().add_edges(vertices + maxCurvDir * length, vertices - maxCurvDir * length, red);
            return true;
        case '3':
            clear_lines(viewer);
            // Draw both segments
            viewer.data().add_edges(vertices + minCurvDir * length, vertices - minCurvDir * length, blue);
            viewer.data().add_edges(vertices + maxCurvDir * length, vertices - maxCurvDir * length, red);
            return true;
        case '4':
            // Remove the segments
            clear_lines(viewer);
            return true;
        default:
            break;
    }
    return false;
}

int main(int argc, char *argv[])
{
    std::string filename = "mesh/bunny.off";

    igl::read_triangle_mesh(filename, vertices, faces);

    // Compute curvature directions via quadric fitting
    igl::principal_curvature(vertices, faces, minCurvDir, maxCurvDir, minCurvVal, maxCurvVal);

    // Mean curvature
    Eigen::VectorXd meanCurv = 0.5 * (minCurvVal + maxCurvVal);

    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(vertices, faces);

    // Compute pseudocolor
    Eigen::MatrixXd color;
    igl::parula(meanCurv, true, color);
    viewer.data().set_colors(color);

    // Compute the one-ring neighbors for all the vertices
    //igl::adjacency_list(faces, adjacencyList);

    // Average edge length for sizing
    length = igl::avg_edge_length(vertices, faces);

    // Hide wireframe
    viewer.data().show_lines = false;

    viewer.callback_key_down = &key_down;

    std::cout << "Press '1' for minimal curvature directions." << std:: endl
              << "Press '2' for maximal curvature directions." << std:: endl
              << "Press '2' for minimal and maximal curvature directions." << std:: endl
              << "Press '2' for no curvature direction." << std:: endl;

    viewer.launch();
}
