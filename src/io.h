//
// Created by victorvianna on 12/12/18.
//

#include<iostream>
#include<vector>

#ifndef SYMMETRY_DETECTION_IO_H
#define SYMMETRY_DETECTION_IO_H

template<typename Iterator>
void write_to_file(Iterator beg, Iterator end, std::string filename) {
    std::ofstream fs(filename);
    for (auto it = beg; it != end; it++)
        fs << *it << std::endl;
    fs.close();
}

template<typename T>
std::ostream &operator<<(std::ostream &os, const std::vector <T>& v) {
    for (auto &x : v)
        os << x << " ";
    return os;
}

#endif //SYMMETRY_DETECTION_IO_H
