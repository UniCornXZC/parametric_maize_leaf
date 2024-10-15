#include "leaf_generation.h"

int main() {

    Leaf_Param leaf_param;
    leaf_param.m1 = 2;
    leaf_param.m2 = 3;
    leaf_param.max_width_posi = 0.65; // MWP
    leaf_param.leaf_length = 1.0; // Unit: m
    leaf_param.leaf_width = 0.08; // Unit: m
    leaf_param.leaf_thickness = 0.0002;
    leaf_param.curve_num = 10;

    bool specific_curves = false;
    if(specific_curves){
        // Create a leaf mesh with specific blade undulations
        leaf_param.argv_y = {13.95, -38.7, 40.5, -44.55, 22.95, -18.9, 30.6, -20.7, 37.8, -18.45}; // Unit: degree
        leaf_param.argv_z = {7.2, -4.8, 20.7, -29.4, 28.5, -21, 3.6, -27, 5.4, -0.6, 16.5};
    }else{
        // Create a leaf mesh with random blade undulations with max angle restrictions
        leaf_param.max_angle_y = 45; // Unit: degree
        leaf_param.max_angle_z = 30;
    }

    // Leaf position coordinates in 3D space
    leaf_param.px = 0, leaf_param.py = 0, leaf_param.pz = 0;
    // Leaf rotation angles in 3D space
    leaf_param.rx = 0, leaf_param.ry = 0, leaf_param.rz = 0;
    // If single layer, the output is just a 3D surface mesh. Otherwise, the output is a manifold 3D object.
    leaf_param.is_single_layer = false;

    // Replace folder path with yours
    std::string dest_folder = "C:\\Users\\Xiang\\Desktop\\";
    Leaf leaf = gen_leaf_om(leaf_param, dest_folder);

    float surface_area;

    float area = CGAL::to_double(CGAL::Polygon_mesh_processing::area(leaf.get_mesh()));
    if(leaf_param.is_single_layer){
        surface_area = area * 10000; // Convert square meter to square centimeter

    }else{
        // We only count single side leaf surface area
        surface_area = area * 5000;

    }
    std::cout << "Surface Area: " << surface_area << std::endl;

    return 0;
}
