//
// Created by gabriel on 13/12/18.
//

#include "patch.h"

PatchPair::PatchPair(Transformation _transform, std::vector<std::vector<int>> &oneRing, double _threshold) : transform(_transform) {
    origin.push_back(transform.getOrigin());
    image.push_back(transform.getImage());
    threshold = _threshold;

    for(int i : oneRing[transform.getImage()])
        imageFrontier.insert(i);
}

int PatchPair::closestFrontier(Eigen::MatrixXd point, Eigen::MatrixXd &V) {
    double minDist = 0;
    int minIndex = -1;
    for(std::unordered_set<int>::iterator it = imageFrontier.begin(); it != imageFrontier.end(); it++){
        double dist = (point - V.row(*it)).norm();
        if(it == imageFrontier.begin() || dist < minDist){
            minDist = dist;
            minIndex = *it;
        }
    }
    return minIndex;
}

bool PatchPair::insert(int newVertex, Eigen::MatrixXd &V, std::vector<std::vector<int>> &oneRing) {
    /// change newVertex to its transformation
    Eigen::MatrixXd point = V.row(newVertex);
    Eigen::MatrixXd transformed = transform.apply(point);

    int closest = closestFrontier(transformed, V);

    if(closest != -1 && (transformed - V.row(closest)).norm() < threshold){
        origin.push_back(newVertex);
        image.push_back(closest);

        imageFrontier.erase(closest);
        for(int i : oneRing[closest])
            imageFrontier.insert(i);

        return true;
    }

    return false;
}

int PatchPair::size(){
    return origin.size();
}
