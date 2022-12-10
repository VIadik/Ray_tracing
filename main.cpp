#include <iostream>
#include <cmath>
#include <cassert>
#include <vector>
#include <ctime>

using namespace std;

const int BYTES_PER_PIXEL = 3; /// red, green, & blue
const int FILE_HEADER_SIZE = 14;
const int INFO_HEADER_SIZE = 40;

unsigned char *createBitmapFileHeader(int height, int stride) {
    int fileSize = FILE_HEADER_SIZE + INFO_HEADER_SIZE + (stride * height);

    static unsigned char fileHeader[] = {
            0, 0,     /// signature
            0, 0, 0, 0, /// image file size in bytes
            0, 0, 0, 0, /// reserved
            0, 0, 0, 0, /// start of pixel array
    };

    fileHeader[0] = (unsigned char) ('B');
    fileHeader[1] = (unsigned char) ('M');
    fileHeader[2] = (unsigned char) (fileSize);
    fileHeader[3] = (unsigned char) (fileSize >> 8);
    fileHeader[4] = (unsigned char) (fileSize >> 16);
    fileHeader[5] = (unsigned char) (fileSize >> 24);
    fileHeader[10] = (unsigned char) (FILE_HEADER_SIZE + INFO_HEADER_SIZE);

    return fileHeader;
}

unsigned char *createBitmapInfoHeader(int height, int width) {
    static unsigned char infoHeader[] = {
            0, 0, 0, 0, /// header size
            0, 0, 0, 0, /// image width
            0, 0, 0, 0, /// image height
            0, 0,     /// number of color planes
            0, 0,     /// bits per pixel
            0, 0, 0, 0, /// compression
            0, 0, 0, 0, /// image size
            0, 0, 0, 0, /// horizontal resolution
            0, 0, 0, 0, /// vertical resolution
            0, 0, 0, 0, /// colors in color table
            0, 0, 0, 0, /// important color count
    };

    infoHeader[0] = (unsigned char) (INFO_HEADER_SIZE);
    infoHeader[4] = (unsigned char) (width);
    infoHeader[5] = (unsigned char) (width >> 8);
    infoHeader[6] = (unsigned char) (width >> 16);
    infoHeader[7] = (unsigned char) (width >> 24);
    infoHeader[8] = (unsigned char) (height);
    infoHeader[9] = (unsigned char) (height >> 8);
    infoHeader[10] = (unsigned char) (height >> 16);
    infoHeader[11] = (unsigned char) (height >> 24);
    infoHeader[12] = (unsigned char) (1);
    infoHeader[14] = (unsigned char) (BYTES_PER_PIXEL * 8);

    return infoHeader;
}

void generateBitmapImage(unsigned char *image, int height, int width, const char *imageFileName) {
    int widthInBytes = width * BYTES_PER_PIXEL;

    unsigned char padding[3] = {0, 0, 0};
    int paddingSize = (4 - (widthInBytes) % 4) % 4;
    int stride = (widthInBytes) + paddingSize;

    FILE *imageFile = fopen(imageFileName, "wb");

    unsigned char *fileHeader = createBitmapFileHeader(height, stride);
    fwrite(fileHeader, 1, FILE_HEADER_SIZE, imageFile);

    unsigned char *infoHeader = createBitmapInfoHeader(height, width);
    fwrite(infoHeader, 1, INFO_HEADER_SIZE, imageFile);

    int i;
    for (i = 0; i < height; i++) {
        fwrite(image + ((height - i - 1) * widthInBytes), BYTES_PER_PIXEL, width, imageFile);
        fwrite(padding, 1, paddingSize, imageFile);
    }

    fclose(imageFile);
}

#define double long double

const double INF = 1e9;

class Vector {
public:
    double x, y, z;

    Vector(Vector a, Vector b) : x(b.x - a.x), y(b.y - a.y), z(b.z - a.z) {}

    Vector(double x, double y, double z) : x(x), y(y), z(z) {}

    Vector() = default;

    double length() const {
        return sqrt(x * x + y * y + z * z);
    }

    void normalize(double new_len) {
        double len = length();
        x /= len;
        y /= len;
        z /= len;
        x *= new_len;
        y *= new_len;
        z *= new_len;
    }

    Vector operator+(Vector ot) const {
        return {x + ot.x, y + ot.y, z + ot.z};
    }

    Vector operator-(Vector ot) const {
        return {x - ot.x, y - ot.y, z - ot.z};
    }

    Vector operator*(double k) {
        return {x * k, y * k, z * k};
    }

    bool operator==(Vector ot) {
        return x == ot.x && y == ot.y && z == ot.z;
    }
};

const Vector ZERO = {-INF, -INF, -INF};


ostream &operator<<(ostream &os, const Vector &p) {
    return os << p.x << ' ' << p.y << ' ' << p.z << '\n';
}

istream &operator>>(istream &in, Vector &p) {
    in >> p.x >> p.y >> p.z;
    return in;
}

double cross_product(Vector a, Vector b) {
    return a.y * b.z - a.z * b.y + a.z * b.x - a.x * b.z + a.x * b.y - a.y * b.x;
}

double dot_product(Vector a, Vector b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

double angle(Vector a, Vector b) {
    return (double) abs(atan2(cross_product(a, b), dot_product(a, b)));
}

double to_degrees(double radian) {
    return radian * (180 / M_PI);
}

double to_radians(double degree) {
    return degree / (180 / M_PI);
}

const double EPS = 1e-12;

double get_min_root(double a, double b, double c) {
    double d = b * b - 4 * a * c;
    if (d < 0) {
        return INF;
    }
    return (-b - sqrt(d)) / (2 * a);;
}

class Screen {
public:
    int h, w;
};

class Color {
public:
    int r, g, b;

    Color(int r, int g, int b) : r(r), g(g), b(b) {}

    Color() = default;

    Color(const string &s) {
        r = stoi(s.substr(1, 2), nullptr, 16);
        g = stoi(s.substr(3, 2), nullptr, 16);
        b = stoi(s.substr(5, 2), nullptr, 16);
    }

    void cut() {
        r = min(r, 255);
        g = min(g, 255);
        b = min(b, 255);
    }

    Color operator+(Color ot) {
        Color c = {r + ot.r, g + ot.g, b + ot.b};
        c.cut();
        return c;
    }

    void operator+=(Color ot) {
        r += ot.r;
        g += ot.g;
        b += ot.b;
        cut();
    }

    Color operator*(double k) const {
        Color c = {(int) (r * k), (int) (g * k), (int) (b * k)};
        c.cut();
        return c;
    }

    void operator*=(double k) {
        r = (int) (r * k);
        g = (int) (g * k);
        b = (int) (b * k);
        cut();
    }

    string hex_string() {
        char hexColor[8];
        snprintf(hexColor, sizeof hexColor, "#%02X%02X%02X", r, g, b);
        return hexColor;
    }
};

void set_pixel(int x, int y, Screen screen, Color color, unsigned char *img) {
    int i = (screen.w * y + x) * 3;
    img[i] = color.b;
    img[i + 1] = color.g;
    img[i + 2] = color.r;
}

class Sphere {
public:
    Vector cent;
    double r;
    Color color;

    Sphere(Vector &cent, double r, Color &color) : cent(cent), r(r), color(color) {}

    Sphere() = default;

    Vector get_cross(Vector v) {
        double a = dot_product(v, v);
        double b = -2 * (dot_product(v, cent));
        double c = dot_product(cent, cent) - r * r;
        double t = get_min_root(a, b, c);
        if (t == INF) {
            return ZERO;
        }
        return v * t;
    }

    Vector get_normal(Vector v) const {
        Vector res = v - cent;
        res.normalize(1);
        return res;
    }
};


pair<int, int> ViewportToScreen(Vector p, Screen s) {
    return {(int) (s.w * (p.x + 0.5)), (int) (s.h * (1 - (p.y + 0.5)))};
}

Vector ScreenToViewport(pair<int, int> cords, Screen s) {
    return {(double) cords.first / s.w - 0.5, (double) -cords.second / s.h + 0.5, 1};
}

Vector PointToViewport(Vector &p) {
    if (p.z < 1) {
        return ZERO;
    }
    Vector res(p.x / p.z, p.y / p.z, 1);
    if (abs(res.x) > 0.5 || abs(res.y) > 0.5) {
        return ZERO;
    } else {
        return res;
    }
}

class Light {
public:
    double value;

    Light() = default;

    explicit Light(double value) : value(value) {}

    virtual double get_light(Sphere s, Vector p) = 0;
};

class Ambient_light : public Light {
public:
    Ambient_light() = default;

    Ambient_light(double value) : Light(value) {}

    double get_light(Sphere s, Vector p) override {
        return value;
    }
};

class Point_light : public Light {
public:
    Vector p;

    Point_light(Vector p, double val) : p(p), Light(val) {}

    Point_light() = default;

    double get_light(Sphere s, Vector v) override {
        Vector norm = s.get_normal(v);
        Vector dir = p - v;
        dir.normalize(1);
        return max((double) 0, dot_product(norm, dir)) * value;
    }

};

class Directional_light : public Light {
public:
    Vector dir;

    Directional_light() = default;

    Directional_light(Vector d, double val) : Light(val) {
        dir = d;
        dir.normalize(1);
    }

    double get_light(Sphere s, Vector p) override {
        Vector norm = s.get_normal(p);
        return max((double) 0, -dot_product(norm, dir)) * value;
    }
};

const int MAXN = 1000;

class Scene {
public:
    vector<Sphere> spheres;
    vector<Light *> lights;

    void add_sphere(Sphere sphere) {
        spheres.push_back(sphere);
    }

    void add_light(Light *light) {
        lights.push_back(light);
    }

    void render(Screen screen) {
#ifndef __APPLE__
        unsigned char img[screen.w * screen.h * 3];
#else
        static unsigned char img[MAXN * MAXN * 3];
#endif
        for (int x = 0; x < screen.w; x++) {
            for (int y = 0; y < screen.h; y++) {
                Vector v = ScreenToViewport({x, y}, screen);
                int ind = -1;
                Vector cross = ZERO;
                for (int j = 0; j < (int) spheres.size(); j++) {
                    Vector res = spheres[j].get_cross(v);
                    if (!(res == ZERO) && (v - res).length() < (cross - res).length()) {
                        ind = j;
                        cross = res;
                    }
                }
                double total_light = 0;
                if (ind != -1 && !(cross == ZERO)) {
                    for (auto &light: lights) {
                        total_light += light->get_light(spheres[ind], cross);
                    }
                    set_pixel(x, y, screen, spheres[ind].color * total_light, img);
                } else {
                    set_pixel(x, y, screen, {0, 0, 0}, img);
                }
            }
        }
        generateBitmapImage(img, screen.h, screen.w, "output.bmp");
    }
};


int main() {
    cout.precision(20);
#ifndef __APPLE__
    Screen screen = {300, 400};
#else
    freopen("input.txt", "r", stdin);
    Screen screen = {MAXN, MAXN};
#endif
    int n;
    cin >> n;
    Sphere s;
    Scene scene;
    for (int i = 0; i < n; i++) {
        cin >> s.cent >> s.r >> s.color.r >> s.color.g >> s.color.b;
        scene.add_sphere(s);
    }
    int m;
    cin >> m;
    for (int i = 0; i < m; i++) {
        string type;
        cin >> type;
        if (type == "d") {
            auto *light = new Directional_light;
            cin >> light->dir >> light->value;
            light->dir.normalize(1);
            scene.add_light(light);
        } else {
            auto *light = new Point_light;
            cin >> light->p >> light->value;
            scene.add_light(light);
        }
    }
    scene.add_light(new Ambient_light(0.2));
    auto start = clock();
    scene.render(screen);
    cerr << "time: " << (clock() - start) * 1.0 / CLOCKS_PER_SEC << '\n';
    return 0;
}
/*
1
0 0 5 1 255 0 0
1
d 1 -1 0 1
*/
