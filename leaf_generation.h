//
// Created by Zhaocheng Xiang (zxiang2@unl.edu) on 11/6/2023.
//
#include <iostream>
#include <cmath>
#include <utility>
#include <vector>
#include <Eigen/Dense>
#include <gsl/gsl_spline.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>


typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef Mesh::Vertex_index vertex_descriptor;
typedef CGAL::Aff_transformation_3<K> aff_transform;

using Eigen::MatrixXd;

#ifndef CPP_CGAL_LEAF_GENERATION_H
#define CPP_CGAL_LEAF_GENERATION_H

#define PI 3.1415926
#define SEG_NUM 200
#define SINGLE_LAYER_FRAME_SIZE 200

struct Leaf_Param {
    float m1, m2;
    float leaf_length;
    float leaf_width, max_width_posi;
    float leaf_thickness;
    int curve_num;
    float max_angle_y, max_angle_z; //degree
    std::vector<float> argv_y;
    std::vector<float> argv_z;
    float px, py, pz;
    float rx, ry, rz;
    bool is_single_layer;
};

struct Midrib {
    std::vector<Point_3> points;
    std::vector<float> angleX; // radian angle
    float length;
};

struct CrossSection {
    std::vector<Point_3> points;
    float ratio;
    float angleX; // pitch
    float angleY; // yaw
    float angleZ; // roll
    float max_width;
};

struct Blade {
    std::vector<float> ratio;
    std::vector<float> angleX; // pitch
    std::vector<float> angleY; // yaw
    std::vector<float> angleZ; // roll
    float max_width;
    float max_width_position;
};



class Leaf {
    Midrib midrib;
    Blade blade;
    std::vector<CrossSection> frames;
    Mesh mesh;
    std::string name;
    Leaf_Param leaf_params;
    float aspect_ratio;


public:
    void create_midrib(float, float, float);

public:
    void create_frames(float);

public:
    void create_blade(float, float, float, int, float, float);

public:
    Mesh create_leaf(float, float, float, float, float, float);

public:
    std::vector<Point_3> get_points();

public:
    std::vector<CrossSection> get_frames();

public:
    void set_name(std::string name_);

public:
    std::string get_name();

public:
    Mesh get_mesh();

public:
    void set_param(Leaf_Param);

public:
    Leaf_Param get_param();

public:
    void create_blade_specific(float leaf_length, float leaf_width, float max_width_position, int curve_num,
                               std::vector<float> argv_y, std::vector<float> argv_z);

public:
    float get_aspect_ratio();
};



void Leaf::set_param(Leaf_Param lp) {
    leaf_params = lp;
}

// More segmentation
Mesh Leaf::create_leaf(float x, float y, float z, float r0, float r1, float r2) {
    leaf_params.px = x, leaf_params.py = y, leaf_params.pz = z;
    leaf_params.rx = r0, leaf_params.ry = r1, leaf_params.rz = r2;

    Mesh m;

    // Make the leaf mesh closed by filling the first and last frame holes if the leaf is double layer
    if (!leaf_params.is_single_layer) {
        for (int i = 0; i < SINGLE_LAYER_FRAME_SIZE - 1; i++) {
            for (int j: {0, SEG_NUM - 1}) {
                vertex_descriptor p0 = m.add_vertex(frames[j].points[i]);
                vertex_descriptor p1 = m.add_vertex(frames[j].points[i + 1]);
                vertex_descriptor q0 = m.add_vertex(frames[j].points[SINGLE_LAYER_FRAME_SIZE * 2 - 1 - i]);
                vertex_descriptor q1 = m.add_vertex(frames[j].points[SINGLE_LAYER_FRAME_SIZE * 2 - 2 - i]);
                if (j == 0) {
                    m.add_face(p0, q1, q0);
                    m.add_face(p0, p1, q1);
                } else {
                    m.add_face(p0, q0, q1);
                    m.add_face(p0, q1, p1);
                }

            }

        }
    }


    // Adjust the i boundary to get leaf segment around different position
    for (int i = 0; i < frames.size() - 1; i++) {
        CrossSection inter0 = frames[i];
        CrossSection inter1 = frames[i + 1];
        int step_num;
        if (leaf_params.is_single_layer) {
            step_num = inter0.points.size() / 2 - 1;

            for (int j = 0; j < step_num; j++) {
                vertex_descriptor p0 = m.add_vertex(inter0.points[j]);
                vertex_descriptor p1 = m.add_vertex(inter0.points[j + 1]);

                vertex_descriptor q0 = m.add_vertex(inter1.points[j]);
                vertex_descriptor q1 = m.add_vertex(inter1.points[j + 1]);

                m.add_face(p0, q1, q0);
                m.add_face(p0, p1, q1);
            }

        } else {
            step_num = inter0.points.size() - 1;

            for (int j = 0; j < step_num; j++) {
                vertex_descriptor p0 = m.add_vertex(inter0.points[j]);
                vertex_descriptor p1 = m.add_vertex(inter0.points[j + 1]);

                vertex_descriptor q0 = m.add_vertex(inter1.points[j]);
                vertex_descriptor q1 = m.add_vertex(inter1.points[j + 1]);

                m.add_face(p0, q0, q1);
                m.add_face(p0, q1, p1);

            }

            vertex_descriptor p0 = m.add_vertex(inter0.points[0]);
            vertex_descriptor p1 = m.add_vertex(inter0.points[SINGLE_LAYER_FRAME_SIZE * 2 - 1]);

            vertex_descriptor q0 = m.add_vertex(inter1.points[0]);
            vertex_descriptor q1 = m.add_vertex(inter1.points[SINGLE_LAYER_FRAME_SIZE * 2 - 1]);

            m.add_face(p0, q1, q0);
            m.add_face(p0, p1, q1);
        }


    }

    glm::mat4 mat = glm::mat4(1.0f);
    glm::mat4 rot_x = glm::rotate(mat, glm::radians(leaf_params.rx), glm::vec3(1.0, 0.0, 0.0));
    glm::mat4 rot_y = glm::rotate(mat, glm::radians(leaf_params.ry), glm::vec3(0.0, 1.0, 0.0));
    glm::mat4 rot_z = glm::rotate(mat, glm::radians(leaf_params.rz), glm::vec3(0.0, 0.0, 1.0));
    glm::mat4 tm = rot_x * rot_y * rot_z;
    aff_transform transformation(tm[0][0], tm[0][1], tm[0][2], leaf_params.px,
                                 tm[1][0], tm[1][1], tm[1][2], leaf_params.py,
                                 tm[2][0], tm[2][1], tm[2][2], leaf_params.pz,
                                 1.0);
    for (Mesh::Vertex_index v: m.vertices()) {
        Point_3 p1 = m.point(v);
        Point_3 p2 = transformation(p1);
        m.point(v) = p2;
    }


    mesh = m;
    return m;
}

// Create a specific length leaf
void Leaf::create_midrib(float m1, float m2, float leaf_length) {

    auto OperatorVec = [](float m1, float m2, float i, float step) {
        float j = i + step;
        float x1 = m1 * cos(i) / i;
        float x2 = m1 * cos(j) / j;
        float y1 = -abs(m2) * sin(i) / i;
        float y2 = -abs(m2) * sin(j) / j;

        double dx = (x2 - x1) / step;
        double dy = (y2 - y1) / step;


        std::vector<double> tmpVec = {x1, y1, dy / dx, sqrt(pow(dx, 2) + pow(dy, 2))};
        return tmpVec;
    };

    float low_bound = 3.0;
    float high_bound = 5.5;
    float step = (high_bound - low_bound) / SEG_NUM;
    std::vector<double> opStart = OperatorVec(m1, m2, low_bound, step);
    std::vector<double> opEnd = OperatorVec(m1, m2, high_bound, step);

    float rib_length = opStart[3] + opEnd[3];

    // Using Trapezoidal Rule to calculate midrib length
    for (int i = 1; i < SEG_NUM; i++) {
        std::vector<double> op = OperatorVec(m1, m2, low_bound + i * step, step);
        rib_length += 2 * op[3];
    }
    rib_length *= step / 2;

    // Scale midrib by defined leaf length
    float scale_ratio = leaf_length / rib_length;
    float y_offset = 0 - opStart[0] * scale_ratio;
    float z_offset = 0 - opStart[1] * scale_ratio;
    Point_3 point = {0, 0, 0};
    midrib.points.push_back(point);
    float angle = atan(opStart[2]) * 180.0f / M_PI;
    midrib.angleX.push_back(angle);
    float rib_length_scaled = (opStart[3] + opEnd[3]) * scale_ratio;
    float highest = 0;
    for (int i = 1; i < SEG_NUM; i++) {
        std::vector<double> op = OperatorVec(m1, m2, low_bound + i * step, step);
        float y = y_offset + op[0] * scale_ratio;
        float z = z_offset + op[1] * scale_ratio;
        if (z > highest)
            highest = z;
        Point_3 temp_point = {0, y, z};
        midrib.points.push_back(temp_point);
        float tmpAngle = atan(op[2]) * 180.0f / M_PI;
        midrib.angleX.push_back(tmpAngle);
        rib_length_scaled += 2 * op[3] * scale_ratio;
    }
    rib_length_scaled *= step / 2;
    midrib.length = rib_length_scaled;
    float y_distance = CGAL::to_double(midrib.points[SEG_NUM - 1].y() - midrib.points[0].y());
    float z_distance = highest - CGAL::to_double(midrib.points[0].z());
    aspect_ratio = z_distance / y_distance;
}

void Leaf::create_frames(float leaf_thickness) {
    if (blade.ratio.empty())
        std::cerr << "Leaf Generation: Ratio vector in Blade is empty!" << std::endl;
    if (blade.ratio.size() != midrib.points.size())
        std::cerr << "Leaf Generation: Ratio vector dimension isn't equal to Blade points dimension!" << std::endl;
    int segment_num = midrib.points.size();
    float low_bound = -5;
    float high_bound = 5;
    int num_points = SINGLE_LAYER_FRAME_SIZE;
    float step = (high_bound - low_bound) / num_points;

    for (int i = 0; i < segment_num; i++) {
        CrossSection frame;
        frame.max_width = blade.max_width;
        frame.ratio = blade.ratio[i] * blade.max_width;
        std::vector<Point_3> upper_bound;
        Point_3 midrib_pt = midrib.points[i];

        // Rotate frame to get blade curves
        float angle_x = midrib.angleX[i] * 1.0f;
        float angle_y = blade.angleY[i] * 1.0f;
        float angle_z = blade.angleZ[i] * 1.0f;

        glm::mat4 trans = glm::mat4(1.0f);
        glm::mat4 rx = glm::rotate(trans, glm::radians(angle_x), glm::vec3(1.0, 0.0, 0.0));
        glm::mat4 ry = glm::rotate(trans, glm::radians(angle_y), glm::vec3(0.0, 1.0, 0.0));
        glm::mat4 rz = glm::rotate(trans, glm::radians(angle_z), glm::vec3(0.0, 0.0, 1.0));
        glm::mat4 rot_mat = rx * ry * rz;

        for (int j = 0; j < num_points; j++) {
            float t = low_bound + j * step;
            float tmp_x = blade.ratio[i] * blade.max_width * (t - tanh(t)) / 8;
            float tmp_z = blade.ratio[i] * blade.max_width * (1 - 1 / cosh(t)) / 4;

            glm::vec4 low_point_vec = {tmp_x, 0, tmp_z, 1};
            glm::vec4 high_point_vec = {tmp_x, 0, tmp_z + leaf_thickness, 1};

            glm::vec4 lp_vec = rot_mat * low_point_vec;
            glm::vec4 hp_vec = rot_mat * high_point_vec;
            Point_3 low_point = {midrib_pt.x() + lp_vec[0], midrib_pt.y() + lp_vec[1], midrib_pt.z() + lp_vec[2]};
            Point_3 high_point = {midrib_pt.x() + hp_vec[0], midrib_pt.y() + hp_vec[1], midrib_pt.z() + hp_vec[2]};

            // Move upward to get leaf thickness
            frame.points.push_back(low_point);
            upper_bound.push_back(high_point);
            // TODO Non-uniform leaf thickness fitting
        }
        for (int k = 1; k <= num_points; k++) {
            frame.points.push_back(upper_bound[num_points - k]);
        }


        frames.push_back(frame);

    }

}

void Leaf::create_blade_specific(float leaf_length, float leaf_width, float max_width_position, int curve_num,
                                 std::vector<float> argv_y, std::vector<float> argv_z) {
    blade.max_width = leaf_width;
    blade.max_width_position = max_width_position;
    float shape_constant = 0.6;

    float low_bound = 0;
    float high_bound = leaf_length;
    float step = (high_bound - low_bound) / SEG_NUM;
    float sum = 0.001;
    blade.ratio.push_back(0.001);
    for (int i = 1; i < SEG_NUM; i++) {
        float width_ratio = pow(sin(((PI * (low_bound + step * i)) / (2 * max_width_position * leaf_length))),
                                shape_constant);
        sum += width_ratio;
        blade.ratio.insert(blade.ratio.begin(), width_ratio);

    }

    // Add X Y Z angle to intersection
    double yi, zi;
    auto *x = new double[curve_num + 1];
    auto *y = new double[curve_num + 1];
    auto *z = new double[curve_num + 1];
    x[curve_num] = SEG_NUM, y[curve_num] = 0, z[curve_num] = 0;
    int interval = SEG_NUM / curve_num;
    std::vector<int> segments;
    std::vector<float> amplitude;
    int cnt = 0;
    for (int i = 0; i < curve_num; i++) {
        float tmpSum = 0;

        for (int j = 0; j < interval; j++) {
            tmpSum += blade.ratio[i * interval + j];
        }
        x[i] = cnt;
        y[i] = argv_y[i];
        z[i] = argv_z[i];

        amplitude.push_back(tmpSum / interval);
        cnt += tmpSum / sum * SEG_NUM;
        segments.push_back(tmpSum / sum * SEG_NUM);
    }


    // Curve Fitting for Blades
    {
        gsl_interp_accel *acc
                = gsl_interp_accel_alloc();

        gsl_spline *spline
                = gsl_spline_alloc(gsl_interp_cspline, curve_num + 1);

        gsl_spline_init(spline, x, y, curve_num + 1);

        for (int a = 0; a < SEG_NUM; a++) {
            yi = gsl_spline_eval(spline, a, acc);
            blade.angleY.push_back(yi);
        }

        gsl_spline_init(spline, x, z, curve_num + 1);

        for (int b = 0; b < SEG_NUM; b++) {
            zi = gsl_spline_eval(spline, b, acc);
            blade.angleZ.push_back(zi);
        }

        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);
    }
}

void Leaf::create_blade(float leaf_length, float leaf_width, float max_width_position, int curve_num,
                        float max_curve_angle_y, float max_curve_angle_z) {
    blade.max_width = leaf_width;
    blade.max_width_position = max_width_position;
    float shape_constant = 0.6;

    float low_bound = 0;
    float high_bound = leaf_length;
    float step = (high_bound - low_bound) / SEG_NUM;
    float sum = 0.001;
    blade.ratio.push_back(0.001);
    for (int i = 1; i < SEG_NUM; i++) {
        float width_ratio = pow(sin(((PI * (low_bound + step * i)) / (2 * max_width_position * leaf_length))),
                                shape_constant);
        sum += width_ratio;
        blade.ratio.insert(blade.ratio.begin(), width_ratio);

    }

    // Add X Y Z angle to intersection
    double yi, zi;
    auto *x = new double[curve_num + 1];
    auto *y = new double[curve_num + 1];
    auto *z = new double[curve_num + 1];
    x[0] = 0;
    x[curve_num] = SEG_NUM, y[curve_num] = 0, z[curve_num] = 0;

    int interval = SEG_NUM / curve_num;
    std::vector<int> segments;
    std::vector<float> amplitude;
    int cnt = 0;
    for (int i = 0; i < curve_num; i++) {
        float tmpSum = 0;

        for (int j = 0; j < interval; j++) {
            tmpSum += blade.ratio[i * interval + j];
        }
        x[i] = cnt;
//        y[i] = tmpSum/interval;
        // add randomness to leaf blade curve
        int flag = 1;
        if (i % 2 == 1) {
            flag = -1;
        }

        y[i] = flag * max_curve_angle_y * (rand() % 100) / 100;
        z[i] = flag * max_curve_angle_z * (rand() % 100) / 100;
        leaf_params.argv_y.push_back(y[i]);
        leaf_params.argv_z.push_back(z[i]);

        amplitude.push_back(tmpSum / interval);
        cnt += tmpSum / sum * SEG_NUM;
        segments.push_back(tmpSum / sum * SEG_NUM);
    }

    // Curve Fitting for Blades
    {
        gsl_interp_accel *acc
                = gsl_interp_accel_alloc();

        gsl_spline *spline
                = gsl_spline_alloc(gsl_interp_cspline, curve_num + 1);

        gsl_spline_init(spline, x, y, curve_num + 1);

        for (int a = 0; a < SEG_NUM; a++) {
            yi = gsl_spline_eval(spline, a, acc);
            blade.angleY.push_back(yi);
        }

        gsl_spline_init(spline, x, z, curve_num + 1);

        for (int b = 0; b < SEG_NUM; b++) {
            zi = gsl_spline_eval(spline, b, acc);
            blade.angleZ.push_back(zi);
        }

        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);
    }


}

std::vector<CrossSection> Leaf::get_frames() { return frames; }

void Leaf::set_name(std::string name_) {
    name = name_;
}

std::string Leaf::get_name() {
    return name;
}

Mesh Leaf::get_mesh() {
    return mesh;
}

float Leaf::get_aspect_ratio() {
    return aspect_ratio;
}

Leaf_Param Leaf::get_param() {
    return leaf_params;
}

std::string gen_random(const int len) {
    static const char alphanum[] =
            "0123456789"
            "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    std::string tmp_s;
    tmp_s.reserve(len);

    for (int i = 0; i < len; ++i) {
        tmp_s += alphanum[rand() % (sizeof(alphanum) - 1)];
    }

    return tmp_s;
}


Leaf gen_leaf_om(Leaf_Param leaf_params, std::string folder) {

    std::string name =
            "leaf_" + std::format("{:.1f}", leaf_params.m1) + "_" + std::format("{:.1f}", leaf_params.m2) + "_" +
            std::format("{:.1f}", leaf_params.leaf_length) + "_" +
            std::format("{:.2f}", leaf_params.leaf_width) + "_" +
            std::format("{:.2f}", leaf_params.max_width_posi) + "-" +
            std::to_string(leaf_params.curve_num) + "_" + std::format("{:.1f}", leaf_params.max_angle_y) + "_" +
            std::format("{:.1f}", leaf_params.max_angle_z) + "-" + std::format("{:.1f}", leaf_params.pz) + "_" +
            std::format("{:.1f}", leaf_params.rx) + "_" +
            std::format("{:.1f}", leaf_params.ry) + "_" + std::format("{:.1f}", leaf_params.rz) + "_" + gen_random(5);
    std::string output_file = name + ".stl";
    Leaf leaf = Leaf();
    leaf.set_param(leaf_params);
    leaf.set_name(name);
    leaf.create_midrib(leaf_params.m1, leaf_params.m2, leaf_params.leaf_length);
    if (leaf_params.argv_y.size() >= leaf_params.curve_num && leaf_params.argv_z.size() >= leaf_params.curve_num) {
        leaf.create_blade_specific(leaf_params.leaf_length, leaf_params.leaf_width, leaf_params.max_width_posi,
                                   leaf_params.curve_num, leaf_params.argv_y, leaf_params.argv_z);
    } else {
        leaf.create_blade(leaf_params.leaf_length, leaf_params.leaf_width, leaf_params.max_width_posi,
                          leaf_params.curve_num, leaf_params.max_angle_y,
                          leaf_params.max_angle_z);
    }

    leaf.create_frames(leaf_params.leaf_thickness);
    Mesh m = leaf.create_leaf(leaf_params.px, leaf_params.py, leaf_params.pz, leaf_params.rx, leaf_params.ry,
                              leaf_params.rz);

    if (!folder.empty()) {
        CGAL::IO::write_polygon_mesh(folder + output_file, m, CGAL::parameters::stream_precision(17));
    }

    return leaf;
}

#endif //CPP_CGAL_LEAF_GENERATION_H
