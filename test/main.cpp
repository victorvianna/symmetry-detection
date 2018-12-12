#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/adjacency_list.h>
#include <igl/avg_edge_length.h>
#include <igl/parula.h>
#include <igl/principal_curvature.h>
#include <igl/read_triangle_mesh.h>
#include "mean_shift.h"
#include "kernel.h"
#include "flann/flann.hpp"

Eigen::MatrixXd vertices;
Eigen::MatrixXi faces;

Eigen::MatrixXd minCurvDir, maxCurvDir, normal;
Eigen::VectorXd minCurvVal, maxCurvVal;

const Eigen::RowVector3d red(0.8, 0.2, 0.2), green(0.2, 0.8, 0.2), blue(0.2, 0.2, 0.8);
double length;

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

// Remove the segment lines added to the mesh
void clear_lines(igl::opengl::glfw::Viewer& viewer){
    viewer.data().lines.resize(0, 9);
}

// This function is called every time a keyboard button is pressed
bool key_down_curvatures(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier){
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
            clear_lines(viewer);
            // Draw normals
            viewer.data().add_edges(vertices + normal * length, vertices - normal * length, green);
            return true;
        case '5':
            clear_lines(viewer);
            // Draw the three directions
            viewer.data().add_edges(vertices + minCurvDir * length, vertices - minCurvDir * length, blue);
            viewer.data().add_edges(vertices + maxCurvDir * length, vertices - maxCurvDir * length, red);
            viewer.data().add_edges(vertices + normal * length, vertices - normal * length, green);
            return true;
        case '6':
            // Remove the segments
            clear_lines(viewer);
            return true;
        default:
            break;
    }
    return false;
}

void test_principal_curvatures(){
    std::string filename = "mesh/bunny.off";

    igl::read_triangle_mesh(filename, vertices, faces);

    // Compute curvature directions via quadric fitting
    igl::principal_curvature(vertices, faces, minCurvDir, maxCurvDir, minCurvVal, maxCurvVal);

    // Compute the normal directions
    /*normal = Eigen::MatrixXd(vertices.rows(), 3);
    for (int i = 0; i < vertices.rows(); i++) {
        Eigen::Vector3d d1, d2, n;
        d1 << minCurvDir.row(i);
        d2 << maxCurvDir.row(i);
        n = d1.cross(d2);
        n.normalize();
        normal.row(i) << n(0), n(1), n(2);
    }*/
    igl::per_vertex_normals(vertices, faces, normal);

    // Mean curvature
    Eigen::VectorXd meanCurv = 0.5 * (minCurvVal + maxCurvVal);

    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(vertices, faces);

    // Compute pseudocolor
    Eigen::MatrixXd color;
    igl::parula(meanCurv, true, color);
    viewer.data().set_colors(color);

    // Average edge length for sizing
    length = igl::avg_edge_length(vertices, faces);

    // Hide wireframe
    viewer.data().show_lines = false;

    viewer.callback_key_down = &key_down_curvatures;

    std::cout << "Press '1' for minimal curvature directions." << std:: endl
              << "Press '2' for maximal curvature directions." << std:: endl
              << "Press '3' for minimal and maximal curvature directions." << std:: endl
              << "Press '4' for normal directions." << std:: endl
              << "Press '5' for the three directions." << std:: endl
              << "Press '6' for no curvature direction." << std:: endl;

    viewer.launch();
}

void test_clustering(){
    using namespace std;
    MeanShift *msp = new MeanShift(epanechnikov_kernel);
    double kernel_bandwidth = 3;

    vector<vector<double> > points = {{0, 0, 0}, {1, 1, 1}};
    vector<Cluster> clusters = msp->cluster(points, kernel_bandwidth);
}

void test_weighted_clustering(){
    using namespace std;
    MeanShift *msp = new MeanShift(epanechnikov_kernel, {1, 1, 1, 0, 0});
    double kernel_bandwidth = 3;

    vector<vector<double> > points = {{0, 0, 0, 10, 20}, {1, 1, 1, 30, 40}};
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
    double* query_ = new double[4]{
            0.5, 0.5,
            2.5, 2.5
    };
    Matrix<double> query(query_, 2, 2);

    std::vector< std::vector<int>> indices;
    std::vector<std::vector<double>> dists;

    // construct an randomized kd-tree index using 4 kd-trees

    Index<L2<double>> index(dataset, flann::KDTreeIndexParams(4));
    index.buildIndex();

    // search around a certain radius
    double radius = 0.6;
    index.radiusSearch(query, indices, dists, radius, flann::SearchParams(128));

    for(std::vector<int>& v : indices) {
        for (int x : v)
            std::cout << x << " ";
        std::cout<<std::endl;
    }
    delete[] dataset.ptr();
    delete[] query.ptr();
}

int main(int argc, char *argv[])
{
    //test_mesh_display();
    test_principal_curvatures();
    test_clustering();
    test_nearest_neighbors();
    test_weighted_clustering();
    return 0;
}
