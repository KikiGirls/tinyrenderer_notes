#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "geometry.h"
#include "model.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);
Model *model = nullptr;
const int width = 800;
const int height = 800;

void line(Vec2i p0, Vec2i p1, TGAImage &image, TGAColor color) {
    bool steep = false;
    if (std::abs(p0.x - p1.x) < std::abs(p0.y - p1.y)) {
        std::swap(p0.x, p0.y);
        std::swap(p1.x, p1.y);
        steep = true;
    }
    if (p0.x > p1.x) {
        std::swap(p0, p1);
    }

    for (int x = p0.x; x <= p1.x; x++) {
        float t = (x - p0.x) / (float) (p1.x - p0.x);
        int y = p0.y * (1. - t) + p1.y * t;
        if (steep) {
            image.set(y, x, color);
        } else {
            image.set(x, y, color);
        }
    }
}

void findbox(Vec3f *t_screen_coords, TGAImage &image, std::vector<Vec2i> &box) {
    // 临时变量
    Vec2i box_min(std::min({t_screen_coords[0].x, t_screen_coords[1].x, t_screen_coords[2].x}),
                  std::min({t_screen_coords[0].y, t_screen_coords[1].y, t_screen_coords[2].y}));
    Vec2i box_max(std::max({t_screen_coords[0].x, t_screen_coords[1].x, t_screen_coords[2].x}),
                  std::max({t_screen_coords[0].y, t_screen_coords[1].y, t_screen_coords[2].y}));

    // 限制在图像范围内
    box_min.x = std::max(0, box_min.x);
    box_min.y = std::max(0, box_min.y);
    box_max.x = std::min(image.get_width() - 1, box_max.x);
    box_max.y = std::min(image.get_height() - 1, box_max.y);

    box = {box_min, box_max};
}


std::array<Vec2f, 3> xyz_to_xy(Vec3f *vec) {
    return {Vec2f(vec[0].x, vec[0].y), Vec2f(vec[1].x, vec[1].y), Vec2f(vec[2].x, vec[2].y)};
}

void getBarycentricIndex(const std::array<Vec2f, 3> &vecs, Vec2i p, Vec3f &bc_index) {
    // 提取三角形的三个顶点
    Vec2f A = vecs[0];
    Vec2f B = vecs[1];
    Vec2f C = vecs[2];

    // 计算三角形的有向面积（用叉积）
    float S_ABC = (B.x - A.x) * (C.y - A.y) - (C.x - A.x) * (B.y - A.y);

    // 如果面积为 0，说明是退化三角形，无法计算重心坐标
    if (std::abs(S_ABC) < 1e-6) {
        bc_index = Vec3f(-1, -1, -1); // 用特殊值表示无效
        return;
    }

    // 计算各重心坐标对应的部分面积
    float S_PBC = (B.x - p.x) * (C.y - p.y) - (C.x - p.x) * (B.y - p.y); // 点 P 和边 BC 的子三角形面积
    float S_PCA = (C.x - p.x) * (A.y - p.y) - (A.x - p.x) * (C.y - p.y); // 点 P 和边 CA 的子三角形面积

    // 归一化，得到重心坐标
    bc_index.x = S_PBC / S_ABC; // 对应顶点 A
    bc_index.y = S_PCA / S_ABC; // 对应顶点 B
    bc_index.z = 1.0f - bc_index.x - bc_index.y; // 对应顶点 C
}

bool inTriangle(Vec3f *t_screen_coords, Vec2i p, Vec3f &bc_index) {
    auto t_screen_coords_xy = xyz_to_xy(t_screen_coords);
    getBarycentricIndex(t_screen_coords_xy, p, bc_index);
    const float epsilon = 1e-6;
    return bc_index[0] >= -epsilon && bc_index[1] >= -epsilon && bc_index[2] >= -epsilon;
}


 void write_triangle(Vec3f *t_screen_coords, std::vector<float> &zBuffer, TGAImage &image, TGAColor color) {
    // 画出轮廓
    // line(t0, t1, image, color);
    // line(t1, t2, image, color);
    // // line(t2, t0, image, color);
    //
    // // 第一种填充三角形的方法
    // // 使用线扫描填充三角形
    // // 就是从顶部一条一条的扫描下来
    // if (t0.y == t1.y && t1.y == t2.y) {
    //     return;
    // }
    // if (t0.y < t1.y) std::swap(t0, t1);
    // if (t0.y < t2.y) std::swap(t2, t0);
    // if (t2.y > t1.y) std::swap(t2, t1);
    //
    // // 对三角顶点的高度进行排序
    // // 对上半部分进行遍历
    //
    // int total_height = t0.y - t2.y;
    //
    // for (int i = 0; i < total_height; i++) {
    //     bool islow = i > t0.y - t1.y || t0.y == t1.y;
    //     float a = float(i) / float(total_height); //比例系数
    //     Vec2i A = t0 + (t2 - t0) * a;
    //
    //     float b = islow ? float(total_height - i) / (t1.y - t2.y) : float(i) / (t0.y - t1.y);
    //     Vec2i B = islow ? t2 + (t1 - t2) * b : t0 + (t1 - t0) * b;
    //
    //     if (A.x > B.x) std::swap(A, B);
    //     // line(A, B, image, color);
    //     // 这样写会导致精度问题然后会有洞
    //     for (int j = A.x; j <= B.x; j++) {
    //         image.set(j, t0.y - i, color);
    //     }
    // }

    // 第二种填充三角形的方法
    // 重心法

    // std::vector<Vec2i> box;
    // findbox(t0, t1, t2, image, box);
    // for (int x = box[0].x; x <= box[1].x; x++) {
    //     for (int y = box[0].y; y <= box[1].y; y++) {
    //         if (inTriangle(t0, t1, t2, x, y)) {
    //             image.set(x, y, color);
    //         }
    //     }
    // }

    // 使用重心坐标来插值计算zbuffer

    std::vector<Vec2i> box;
    findbox(t_screen_coords, image, box);

    for (int y = box[0].y; y <= box[1].y; y++) {
        int row_offset = y * width; // 提前计算行偏移量
        for (int x = box[0].x; x <= box[1].x; x++) {
            Vec2i p = Vec2i(x, y);
            Vec3f bc_index;

            // 检查点是否在三角形内并计算重心坐标
            if (inTriangle(t_screen_coords, p, bc_index)) {
                // 使用重心坐标插值计算深度值
                float p_z = 0;
                for (int i = 0; i < 3; i++) {
                    p_z += t_screen_coords[i].z * bc_index[i];
                }

                // 深度缓冲区检查并更新
                int z_index = x + row_offset; // 计算索引
                if (p_z > zBuffer[z_index]) {
                    zBuffer[z_index] = p_z; // 更新深度值
                    image.set(x, y, color); // 绘制像素
                }
            }
        }
    }
}


Vec3f getTriangleNormal(Vec3f *triangle) {
    // 获取三角形的三个顶点
    Vec3f v0 = triangle[0];
    Vec3f v1 = triangle[1];
    Vec3f v2 = triangle[2];

    // 计算两个边向量
    Vec3f edge1 = v2 - v0;
    Vec3f edge2 = v1 - v0;

    // 计算法向量：边向量的叉积
    Vec3f normal = edge1^ edge2;

    // 确保法向量归一化（如果需要单位法向量）
    normal.normalize();

    return normal; // 返回右值引用
}

Vec3f MPV (Vec3f v) {
    return Vec3f(int((v.x+1.)*width/2.+.5), int((v.y+1.)*height/2.+.5), v.z);
}

int main(int argc, char **argv) {
    //读取model文件
    if (2 == argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head/african_head.obj");
    }
    TGAImage image(width, height, TGAImage::RGB);
    std::vector<float> zBuffer(width * height, -std::numeric_limits<float>::max());

    Vec3f light_dir = Vec3f(0, 0, -1);

    for (int i = 0; i < model->nfaces(); i++) {
        std::vector<int> face = model->face(i);
        Vec3f t_screen_coords[3];
        Vec3f t_world_coords[3];
        Vec3f t_normal;
        for (int j = 0; j < 3; j++) {
            Vec3f v = model->vert(face[j]);
            t_screen_coords[j] = MPV(v);
            t_world_coords[j] = v;
        }

        t_normal = getTriangleNormal(t_world_coords);

        float t_intensity = t_normal * light_dir;
        if (t_intensity > 0) {
            //如果小于0就看不见，代表光线在三角形的背面。
            TGAColor t_color = TGAColor(t_intensity * 255, t_intensity * 255, t_intensity * 255, 255);
            write_triangle(t_screen_coords, zBuffer, image, t_color);
        }
    }


    image.flip_vertically();
    image.write_tga_file("output.tga");
    delete model;
    return 0;
}
