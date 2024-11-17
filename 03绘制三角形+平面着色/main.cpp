#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "geometry.h"
#include "model.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);
Model *model = NULL;
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

void findbox(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, std::vector<Vec2i> &box) {
    // 临时变量
    Vec2i box_min(std::min({t0.x, t1.x, t2.x}), std::min({t0.y, t1.y, t2.y}));
    Vec2i box_max(std::max({t0.x, t1.x, t2.x}), std::max({t0.y, t1.y, t2.y}));

    // 限制在图像范围内
    box_min.x = std::max(0, box_min.x);
    box_min.y = std::max(0, box_min.y);
    box_max.x = std::min(image.get_width() - 1, box_max.x);
    box_max.y = std::min(image.get_height() - 1, box_max.y);

    box = {box_min, box_max};
}

bool inTriangle(Vec2i A, Vec2i B, Vec2i C, int x, int y) {
    auto cross = [](Vec2i a, Vec2i b) {
        return a.x * b.y - a.y * b.x;
    };
    Vec2i p(x, y);

    Vec2i AB = B - A,BC = C - B, CA = A - C;
    Vec2i A2p = p - A, B2p = p - B, C2p = p - C;

    return (cross(AB, A2p) >= 0 && cross(BC, B2p) >= 0 && cross(CA, C2p) >= 0) ||
           (cross(AB, A2p) <= 0 && cross(BC, B2p) <= 0 && cross(CA, C2p) <= 0);
}


void triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color) {
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

    std::vector<Vec2i> box;
    findbox(t0, t1, t2, image, box);
    for (int x = box[0].x; x <= box[1].x; x++) {
        for (int y = box[0].y; y <= box[1].y; y++) {
            if (inTriangle(t0, t1, t2, x, y)) {
                image.set(x, y, color);
            }
        }
    }
}


Vec3f triangleNormal(Vec3f *triangle) {
    // 获取三角形的三个顶点
    Vec3f v0 = triangle[0];
    Vec3f v1 = triangle[1];
    Vec3f v2 = triangle[2];

    // 计算两个边向量
    Vec3f edge1 = v2 - v0;
    Vec3f edge2 = v1 - v0;

    // 计算法向量：边向量的叉积
    Vec3f normal = edge1 ^ edge2;

    // 确保法向量归一化（如果需要单位法向量）
    normal.normalize();

    return normal; // 返回右值引用
}

int main(int argc, char **argv) {
    //读取model文件
    if (2 == argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }
    TGAImage image(width, height, TGAImage::RGB);

    Vec3f light_dir = Vec3f(0, 0, -1);

    for (int i = 0; i < model->nfaces(); i++) {
        std::vector<int> face = model->face(i);
        Vec2i screen_coords[3];
        Vec3f world_coords[3];
        Vec3f n;
        for (int j = 0; j < 3; j++) {
            Vec3f v = model->vert(face[j]);
            screen_coords[j] = Vec2i((v.x + 1.) * width / 2., (v.y + 1.) * height / 2.);
            world_coords[j] = v;
        }

        n = triangleNormal(world_coords);

        float intensity = n * light_dir;
        if (intensity > 0) {//如果小于0就看不见，代表光线在三角形的背面。
            TGAColor color = TGAColor(intensity * 255, intensity * 255, intensity * 255, 255);
            triangle(screen_coords[0], screen_coords[1], screen_coords[2], image, color);
        }
    }


    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");
    return 0;
}
