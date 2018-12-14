//
// Created by gabriel on 13/12/18.
//

#include "verification.h"
#include "transformation.h"
#include <igl/adjacency_list.h>
#include <set>
#include <unordered_map>
#include <queue>

Verification::Verification(Eigen::MatrixXi& F, double _threshold) {
    threshold = _threshold;
    igl::adjacency_list(F, oneRing);
}

/// TODO
// create class Cluster

void Verification::verifyCluster(Eigen::MatrixXd& V, std::vector<Transformation> &cluster, std::vector<PatchPair> &patches) {
    std::queue<int> extendPatch;
    std::unordered_map<int, bool> visited;

    int patchLength = patches.size();

    // QUESTION
    // Add point to the patch after of before checking its oneRing?

    for(Transformation &transform : cluster) {
        if(!visited[transform.getOrigin()]) {
            visited[transform.getOrigin()] = true;
            extendPatch.push(transform.getOrigin());

            patchLength++;
            patches.push_back(PatchPair(transform, oneRing, threshold));
        }

        while(!extendPatch.empty()) {
            int vertex = extendPatch.front();
            extendPatch.pop();

            for(int i : oneRing[vertex]) {
                if(!visited[i]) {
                    visited[i] = true;

                    if (patches[patchLength - 1].insert(i, V, oneRing))
                        extendPatch.push(i);
                }
            }
        }
    }
}