#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "model.h"

Model::Model(const char *filename) : verts_(), faces_(), textures_() {
    std::ifstream in;
    in.open (filename, std::ifstream::in);
    if (in.fail()) return;
    std::string line;
    while (!in.eof()) {
        std::getline(in, line);
        std::istringstream iss(line.c_str());
        char trash;
        if (!line.compare(0,2, "v ")) {
            iss >> trash;
            Vec3f v;
            for (int i=0;i<3;i++) iss >> v[i];
            verts_.push_back(v);
        } else if (!line.compare(0, 3, "vt ")) {
            iss >> trash >> trash;
            Vec3f vt;
            for (int i=0;i<3;i++) iss >> vt[i];
            textures_.push_back(vt);
        } else if (!line.compare(0,2, "f ")) {
            std::vector<int> f;
            int itrash, idx, texture_idx;
            iss >> trash;
            while (iss >> idx >> trash >> texture_idx >> trash >> itrash) {
                idx--;
                texture_idx--;
                f.push_back(idx);
                f.push_back(texture_idx);
            }
            faces_.push_back(f);
        }
        //std::cerr << "# v#" << verts_.size() << " vt#" << textures_.size() << " f# " << faces_.size() << std::endl;
    }
}

Model::~Model() {
}

int Model::nverts() {
    return (int)verts_.size();
}

int Model::nfaces() {
    return (int)faces_.size();
}

int Model::ntextures() {
    return (int)textures_.size();
}

std::vector<int> Model::face(int idx) {
    return faces_[idx];
}

Vec3f Model::vert(int i) {
    return verts_[i];
}

Vec3f Model::texture(int i) {
    return textures_[i];
}