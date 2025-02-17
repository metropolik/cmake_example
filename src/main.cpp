#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <math.h>
#include <iostream>
#include <unordered_map>
#include <string>

std::unordered_map<std::string, std::vector<float>> sourceImages;

int add(int i, int j) {
    return i + j;
}

bool hasImage(std::string imgName) {
    return !(sourceImages.find(imgName) == sourceImages.end());
}

void addNewImage(std::string imgName, std::vector<float> img) {
    std::cout << "Adding new image " << imgName << " to collection" << std::endl;
    sourceImages[imgName] = img;
    std::cout << "Done adding" << std::endl;
}

void readSrcPixel(std::string textureName, int srcTexW, int srcTexH, double xco, double yco, float & r, float & g, float & b, float & a) {
    if (!hasImage(textureName)) {
        std::cout << "FATAL image not contained in cached images!" << std::endl;
        return;
    }
    std::vector<float> & sourceTexture = sourceImages[textureName];
    int nxco = (int)round(xco);
    int nyco = (int)round(yco);
    if (nxco < 0 || nxco >= srcTexW) {
        nxco = (srcTexW + (nxco%srcTexW))%srcTexW;
    }
    if (nyco < 0 || nyco >= srcTexH) {
        nyco = (srcTexH + (nyco%srcTexH))%srcTexH;
    }
    nyco -= 1;
    int r_idx = (nyco * srcTexW + nxco) * 4;
    r = sourceTexture[r_idx];
    g = sourceTexture[r_idx+1];
    b = sourceTexture[r_idx+2];
    a = sourceTexture[r_idx+3];
}

double dist1(double a, double b) {
    return abs(a-b);
}

double dist2(double ax, double ay, double bx, double by) {
    const double f1 = ax - bx;
    const double f2 = ay - by;
    return sqrt(f1*f1 + f2 * f2);
}

void calcWeights(double ax, double ay, double bx, double by, double cx, double cy, double px, double py, double & u, double & v, double & w, bool nudge = false) {
        if (nudge) {
            //avoid edge cases by nuding the point a tiny amount to the center
            const double centerTriangleX = (ax + bx + cx)/3.0;
            const double centerTriangleY = (ay + by + cy)/3.0;
            const double relox = px - centerTriangleX;
            const double reloy = py - centerTriangleY;
            px = centerTriangleX + relox * 0.99;
            py = centerTriangleY + reloy * 0.99;
        }

        //calc barycentric weights
        const double v0x = bx - ax;
        const double v0y = by - ay;
        const double v1x = cx - ax;
        const double v1y = cy - ay;
        const double v2x = px - ax;
        const double v2y = py - ay;
        const double d00 = v0x * v0x + v0y * v0y;
        const double d01 = v0x * v1x + v0y * v1y;
        const double d11 = v1x * v1x + v1y * v1y;
        const double d20 = v2x * v0x + v2y * v0y;
        const double d21 = v2x * v1x + v2y * v1y;
        const double denom = d00 * d11 - d01 * d01;
        v = (d11 * d20 - d01 * d21) / denom;
        w = (d00 * d21 - d01 * d20) / denom;
        u = 1.0 - v - w;
        const double ud = dist1(u, 0.5);
        const double vd = dist1(v, 0.5);
        const double wd = dist1(w, 0.5);
        if ((ud > 0.499 && ud < 0.501) || (vd > 0.499 && vd < 0.501) || (wd > 0.499 && wd < 0.501)) {
            calcWeights(ax, ay, bx, by, cx, cy, px, py, u, v, w, true);
        }
}

bool isInTriangle(const double u, const double v, const double w) {
    if (u < 0.0 || u > 1.0) { return false;}
    if (v < 0.0 || v > 1.0) { return false;}
    if (w < 0.0 || w > 1.0) { return false;}
    return true;
}
/*bool isInTriangle(double ax, double ay, double bx, double by, double cx, double cy, double px, double py) {
    double u, v, w;
    calcWeights(&u, &v, &w);
    if (u < 0.0 || u > 1.0) { return false;}
    if (v < 0.0 || v > 1.0) { return false;}
    if (w < 0.0 || w > 1.0) { return false;}
    return true;
}*/


std::vector<float> remap(int ux, int uy, int uw, int uh, int timw, int timh, std::string textureName, int srcTexW, int srcTexH, std::vector<double> originalUv, std::vector<double> newUv) {
    //std::cout << "Starting remap for " << uw * uh << " pixels" << std::endl;
    std::vector<float> out;
    out.reserve(uw * uh * 4);
    for (int pixel_x = ux; pixel_x < ux + uw; pixel_x++) {
        //std::cout << "Process " << ((double)pixel_x/((double)(ux + uw)))*100.0 << std::endl;
        for (int pixel_y = uy; pixel_y < uy + uh; pixel_y++) {
            //find triangle in which current pixel is in
            double uvp_x = ((double)pixel_x)/((double)timw);
            double uvp_y = 1.0 - ((double)pixel_y)/((double)timh);
            double u, v, w, newCox, newCoy;
            //try first triangle
            //                   0                   1                   2
            calcWeights(newUv[0], newUv[1], newUv[2], newUv[3], newUv[4], newUv[5], uvp_x, uvp_y, u, v, w);
            if (!isInTriangle(u, v, w)) {
            //                       2                   3                   0
                calcWeights(newUv[4], newUv[5], newUv[6], newUv[7], newUv[0], newUv[1], uvp_x, uvp_y, u, v, w);
                newCox = u * originalUv[4] + v * originalUv[6] + w * originalUv[0];
                newCoy = u * originalUv[5] + v * originalUv[7] + w * originalUv[1];
            } else {
                newCox = u * originalUv[0] + v * originalUv[2] + w * originalUv[4];
                newCoy = u * originalUv[1] + v * originalUv[3] + w * originalUv[5];
            }
            newCox *= srcTexW;
            newCoy *= srcTexH;
            float r, g, b, a;
            readSrcPixel(textureName, srcTexW, srcTexH, newCox, newCoy, r, g, b, a);
            out.push_back(r);
            out.push_back(g);
            out.push_back(b);
            out.push_back(a);
        }
    }
    //std::cout << "Remap finished" << std::endl;
    return out;
}

namespace py = pybind11;

PYBIND11_MODULE(cmake_example, m) {
    m.doc() = R"pbdoc(
		Know it or leave it.
    )pbdoc";

    m.def("remap", &remap, R"pbdoc(Remaps the thing.)pbdoc");
    m.def("hasImage", &hasImage, R"pbdoc(Checks if the image is already cached.)pbdoc");
    m.def("addNewImage", &addNewImage, R"pbdoc(Checks if the image is already cached.)pbdoc");

    m.def("add", &add, R"pbdoc( Add two numbers Some other explanation about the add function.)pbdoc");

    m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc( Subtract two numbers Some other explanation about the subtract function.)pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
