#include "tgaimage.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);

// 第一次尝试
// void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
//     for (float t=0.; t<1.; t+=.1) {
//         //t为步长，每次增加t*(x1-x0)或者t*(y1-y0)
//         int x = x0*(1.-t) + x1*t;
//         int y = y0*(1.-t) + y1*t;
//         image.set(x, y, color);
//     }
// }
// 这样会导致计算增加，比如我们只需要算13到80，一共67次，但是我们使用t过后就会计算一百次
// 于是我们直接用x来循环

// 第二次尝试
// void line(int x0, int y0, int x1, int y1, TGAImage &image, const TGAColor &color) {
//     for (int x = x0; x <= x1; x++) {
//         float t = (x - x0) / (float) (x1 - x0);
//         int y = y0 + t * (y1 - y0);
//         image.set(x, y, color);
//     }
// }

// 还有几个新的问题
// 如果斜率大于1的时候，我们应该用y来遍历
// 另外也有可能直线的斜率是负数，也就是反向的

// 第三次尝试
void line(int x1, int y1, int x2, int y2, TGAImage &image, TGAColor color) {
    if (x1 > x2) {
        //如果方向不对就交换位置
        std::swap(x1, x2);
        std::swap(y1, y2);
    }
    bool steep = (std::abs(y2 - y1) > std::abs(x2 - x1));
    if (steep) //如果斜率大于1；
    {
        std::swap(x1, y1);
        std::swap(x2, y2);
    }

    for (int x = x1; x <= x2; x++) {
        float t = (x - x1) / (float) (x2 - x1);

        int y = y1 + t * (y2 - y1);
        if (steep) {
            image.set(y, x, color);
        } else {
            image.set(x, y, color);
        }

    }
}

int main(int argc, char **argv) {
    TGAImage image(100, 100, TGAImage::RGB);
    line(13, 20, 80, 40, image, white); //线段A
    line(20, 13, 40, 80, image, red); //线段B
    line(80, 40, 13, 20, image, red);
    image.flip_vertically();
    image.write_tga_file("output.tga");
    return 0;
}
