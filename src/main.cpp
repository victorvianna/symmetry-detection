#include <iostream>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/principal_curvature.h>
#include <igl/bounding_box.h>
#include "flann/flann.hpp"
#include "mean_shift.h"
#include "kernel.h"
#include "transformation.h"
#include "io.h"
#include "signature.h"

using namespace std;

vector<Signature> prune_points(vector<Signature> &signatures) {
    vector<Signature> pruned;

    for (Signature &s : signatures) {
        if (s.is_not_umbilical_point())
            pruned.push_back(s);
    }

    return pruned;
}


template<typename T>
vector<T> random_sample(vector<T> &v, int num_samples) {
    std::vector<int> indices(v.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::random_shuffle(indices.begin(), indices.end());
    vector<T> sampled;
    for (int idx : indices)
        sampled.push_back(v[idx]);
    return sampled;
};

void build_pairing_kd_tree(vector<Signature> &signatures, vector<Transformation> &transf, bool rigid,
                   string filename = "final_transf_space.txt") {
    transf.clear();
    std::ofstream fs(filename);
    const int NUM_SAMPLES = std::min<int>((int) signatures.size(), 100);
    vector<Signature> samples = random_sample(signatures, NUM_SAMPLES);
    flann::Matrix<double> datapoints(Signature::flatten(signatures, rigid), signatures.size(), Signature::dimension(rigid));
    flann::Matrix<double> query(Signature::flatten(samples, rigid), samples.size(), Signature::dimension(rigid));
    flann::Index<flann::L2<double>> index(datapoints, flann::KDTreeIndexParams(4));
    index.buildIndex();
    double radius = 0.6;
    std::vector<std::vector<int>> indices;
    std::vector<std::vector<double>> dists;
    index.radiusSearch(query, indices, dists, radius, flann::SearchParams(128));
    for (int i = 0; i < samples.size(); i++) {
        vector<int> &neighbors = indices[i];
        Signature &p_a = samples[i];
        for (int j : neighbors) {
            Signature &p_b = signatures[j];
	    if (p_a.get_point_index() == p_b.get_point_index()) // avoid identity transformation
		    continue;
            auto t = Transformation(p_a, p_b);
            fs << t << std::endl;
            transf.push_back(t);
            // t = Transformation(p_b, p_a);
            // fs << t << std::endl;
            // transf.push_back(t);
        }
    }
    fs.close();
}

void build_pairing_all(vector<Signature> &signatures,  string filename) {
    std::ofstream fs(filename);
    for (auto p = signatures.begin(); p != signatures.end(); p++) {
        for (auto q = p + 1; q != signatures.end(); q++){
            fs << Transformation(*p, *q).to_point() << std::endl;
            fs << Transformation(*q, *p).to_point() << std::endl;
        }
    }
    fs.close();
}

void run_clustering(vector<Transformation> &transf_space, vector<vector<Transformation>> &clusters_transf,
                    double diagonal_length) {
    clusters_transf.clear();
    // setup coefficients according to paper
    double beta_1 = 0.01;
    double beta_2 = 1.0 / (M_PI * M_PI);
    double beta_3 = 4.0 / (diagonal_length * diagonal_length);

    double kernel_bandwidth = 3;
    vector<double> weights = {beta_1, beta_2, beta_2, beta_2, beta_3, beta_3, beta_3, 0, 0};
    MeanShift *msp = new MeanShift(epanechnikov_kernel, weights);

    vector<vector<double> > points;
    Transformation::to_points(transf_space, points);
    vector<Cluster> clusters = msp->cluster(points, kernel_bandwidth);
    for (Cluster &cluster : clusters) {
        vector<Transformation> cluster_transf;
        for (vector<double> &point : cluster.original_points)
            cluster_transf.push_back(Transformation(point));
        clusters_transf.push_back(cluster_transf);
    }

}

int main() {

    bool rigid = true;
    bool plotting = true;

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::read_triangle_mesh("mesh/bunny.off", V, F);

    cout << "Total vertices " << V.size() << endl;

    // Plot the mesh
    //igl::opengl::glfw::Viewer viewer;
    //viewer.data().set_mesh(V, F);

    // signature computation
    vector<Signature> signatures;
    Signature::build_signatures(V, F, signatures);
    if (plotting)
        build_pairing_all(signatures, "complete_transf_space.txt");

    // pruning
    signatures = prune_points(signatures);
    cout << "Remaining vertices after pruning " << signatures.size() << endl;
    if (plotting)
        build_pairing_all(signatures, "pruned_transf_space.txt");


    // Plot after pruning
    //Signature::plot_all_directions(viewer, V, F, signatures);

    // pairing
    vector<Transformation> transf_space;
    build_pairing_kd_tree(signatures, transf_space, rigid);
    cout << "Size of transformation space " << transf_space.size() << endl;

    // calculate diagonal of bounding box
    Eigen::Vector3d m = V.colwise().minCoeff();
    Eigen::Vector3d M = V.colwise().maxCoeff();
    auto bounding_box = (M - m);
    double diagonal_length = bounding_box.norm();

    // clustering
    vector<vector<Transformation>> clusters;
    run_clustering(transf_space, clusters, diagonal_length);

    cout << "Num of clusters " << clusters.size() << endl;

    // run verification and build patches
    /// TODO

    // display patches or write to file
    /// TODO

    //viewer.launch();

    return 0;
}
