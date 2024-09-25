#include <iostream>
#include <fstream>
#include <cmath>
#include <optional>
#include <vector>
#include <algorithm>
#include <cstdint>


// Structure de Vecteur 3D
struct Vec3 {
    float x, y, z; // coordonnées du vecteur

    // Constructeur explicite pour initialiser un Vec3 avec des valeurs
    Vec3(float x_ = 0, float y_ = 0, float z_ = 0) : x(x_), y(y_), z(z_) {}


    // Opérateurs de vecteurs
    Vec3 operator-(const Vec3& other) const {
        return { x - other.x, y - other.y, z - other.z }; // soustraction de vecteurs par un autre vecteur
    }
    Vec3 operator+(const Vec3& other) const {
        return { x + other.x, y + other.y, z + other.z }; // addition de vecteurs par un autre vecteur
    }
    Vec3 operator*(float scalar) const { 
        return { x * scalar, y * scalar, z * scalar }; // multiplication de vecteurs par un float
    }
    Vec3 operator*(const Vec3& other) const {
        return { x * other.x, y * other.y, z * other.z }; // multiplication de vecteurs par un autre vecteur
    }
    // Opérateur d'indexation pour accéder aux éléments x, y, z via un indice
    float& operator[](int i) {
        if (i == 0) return x;
        else if (i == 1) return y;
        else if (i == 2) return z;
        throw std::out_of_range("Index out of range");
    }

    const float& operator[](int i) const {
        if (i == 0) return x;
        else if (i == 1) return y;
        else if (i == 2) return z;
        throw std::out_of_range("Index out of range");
    }
};


// Produit scalaire
float dot(const Vec3& a, const Vec3& b) {
    return ((a.x * b.x) + (a.y * b.y) + (a.z * b.z));
}

// Calcul de la norme d'un vecteur 3D pour avoir sa longueur
float length(const Vec3& a) {
    return std::sqrt(dot(a, a));
}
// Normalisation d'un vecteur 3D
Vec3 normalize(const Vec3& v) {
    float len = length(v);
    return { v.x / len, v.y / len, v.z / len };
}

// Norme au carré
float length_squared(const Vec3& a) {
    return dot(a, a);
}





// structure de Sphère
struct Sphere {
    Vec3 center;
    float radius;
    Vec3 albedo;
};

// structure de Rayon
struct Ray {
    Vec3 origin;
    Vec3 direction;
};

struct Light {
    Vec3 origin;
    Vec3 emission;
};

struct Intersection {
    Vec3 point;
    Vec3 normal;
    float distance;
    Vec3 albedo;
};


// Structure de Scène
struct Scene {
    std::vector<Sphere> spheres;
    std::unique_ptr<BVHNode> bvh_root;

    Scene(std::vector<Sphere>& objects) {
        bvh_root = construct_bvh(objects, 0, objects.size());
    }

    std::optional<Intersection> intersect(const Ray& ray) const {
        return bvh_root->intersect(ray);
    }
};

/////////////////////////////////////////////////
struct AABB {
    Vec3 min;
    Vec3 max;
    

    // Calcul de l'intersection avec un rayon
    bool intersect(const Ray& ray, float& tmin, float& tmax) const {
        for (int i = 0; i < 3; ++i) {
            float invD = 1.0f / ray.direction[i];
            float t0 = (min[i] - ray.origin[i]) * invD;
            float t1 = (max[i] - ray.origin[i]) * invD;
            if (invD < 0.0f) std::swap(t0, t1);
            tmin = std::max(tmin, t0);
            tmax = std::min(tmax, t1);
            if (tmax <= tmin) return false;
        }
        return true;
    }

    // Fusionner deux boîtes englobantes
    static AABB surrounding_box(const AABB& box0, const AABB& box1) {
        Vec3 small(fmin(box0.min.x, box1.min.x), fmin(box0.min.y, box1.min.y), fmin(box0.min.z, box1.min.z));
        Vec3 big(fmax(box0.max.x, box1.max.x), fmax(box0.max.y, box1.max.y), fmax(box0.max.z, box1.max.z));
        return AABB{ small, big };
    }

};
/////////////////////////////////////////////////





// Fonction d'intersection entre un rayon et une sphère
std::optional<Intersection> intersect_sphere(const Ray& ray, const Sphere& sphere) {
    Vec3 oc = ray.origin - sphere.center; // vecteur entre l'origine du rayon et le centre de la sphère

    float a = length_squared(ray.direction); // norme du vecteur direction du rayon
    float b = 2.0f * dot(oc, ray.direction); // produit scalaire entre le vecteur oc et le vecteur direction du rayon
    float c = length_squared(oc) - sphere.radius * sphere.radius; // norme du vecteur oc - rayon de la sphère

    float delta = b * b - 4.0f * a * c;

    if (delta >= 0.0f) {
        // calcul des deux solutions de l'équation du second degré
        float t1 = (-b - std::sqrt(delta)) / (2.0f * a);
        float t2 = (-b + std::sqrt(delta)) / (2.0f * a);
        float t = (t1 >= 0.0f) ? t1 : t2; // on prend la solution positive

        if (t >= 0.0f) {
            Vec3 intersection_point = ray.origin + ray.direction * t;
            Vec3 normal = normalize(intersection_point - sphere.center);
            return Intersection{ intersection_point, normal, t, sphere.albedo };
        }
    }
    return std::nullopt;
}



////////////////////////////////////////////////////////
// Structure BVHNode (nœud de la hiérarchie de volumes englobants)
struct BVHNode {
    std::unique_ptr<BVHNode> left;  // Fils gauche
    std::unique_ptr<BVHNode> right; // Fils droit
    AABB box;  // Boîte englobante
    std::optional<Sphere> object;  // Objet (feuille)
    

    // Constructeur pour les feuilles
    BVHNode(const Sphere& sphere) : object(sphere) {
        box = AABB{ sphere.center - Vec3{sphere.radius, sphere.radius, sphere.radius},
                   sphere.center + Vec3{sphere.radius, sphere.radius, sphere.radius} };
    }

    // Constructeur pour les nœuds internes
    BVHNode(std::unique_ptr<BVHNode> left, std::unique_ptr<BVHNode> right)
        : left(std::move(left)), right(std::move(right)) {
        box = AABB::surrounding_box(this->left->box, this->right->box);
    }

    // Fonction d'intersection avec le BVH
    std::optional<Intersection> intersect(const Ray& ray) {
        float tmin = 0.0f, tmax = std::numeric_limits<float>::infinity();
        if (!box.intersect(ray, tmin, tmax)) return std::nullopt;

        if (object.has_value()) {
            return intersect_sphere(ray, object.value());
        }

        auto hit_left = left ? left->intersect(ray) : std::nullopt;
        auto hit_right = right ? right->intersect(ray) : std::nullopt;

        if (hit_left && hit_right) {
            return (hit_left->distance < hit_right->distance) ? hit_left : hit_right;
        }
        return hit_left ? hit_left : hit_right;
    }

};

// Construction du BVH à partir des sphères
std::unique_ptr<BVHNode> construct_bvh(std::vector<Sphere>& spheres, int start, int end) {
    int n = end - start;

    if (n == 1) {
        return std::make_unique<BVHNode>(spheres[start]);  // Feuille contenant une sphère
    }

    int axis = rand() % 3;  // Choix aléatoire de l'axe de séparation
    std::sort(spheres.begin() + start, spheres.begin() + end, [axis](const Sphere& a, const Sphere& b) {
        return a.center[axis] < b.center[axis];
        });

    int mid = start + n / 2;
    return std::make_unique<BVHNode>(
        construct_bvh(spheres, start, mid),
        construct_bvh(spheres, mid, end)
    );
}
///////////////////////////////////////////////////////



//// Calcul de l'intersection avec toutes les sphères
//std::optional<Intersection> intersect_scene(const Ray& ray, const std::vector<Sphere>& spheres) {
//    std::optional<Intersection> closest_intersection;
//    for (const auto& sphere : spheres) {
//        auto intersection = intersect_sphere(ray, sphere);
//        if (intersection && (!closest_intersection || intersection->distance < closest_intersection->distance)) {
//            closest_intersection = intersection;
//        }
//    }
//    return closest_intersection;
//}

// Calcul de l'intersection avec toutes les sphères
std::optional<Intersection> intersect_scene(const Ray& ray, const Scene scene) {
    std::optional<Intersection> closest_intersection;
    for (const auto& sphere : scene.spheres) {
        auto intersection = intersect_sphere(ray, sphere);
        if (intersection && (!closest_intersection || intersection->distance < closest_intersection->distance)) {
            closest_intersection = intersection;
        }
    }
    return closest_intersection;
}



// Calcul du tonemapping pour ajuster la couleur
Vec3 tonemap(const Vec3& color, float scale) {
    return { std::min(255.0f, color.x * scale), std::min(255.0f, color.y * scale), std::min(255.0f, color.z * scale) };
}

//// Fonction principale pour lancer les rayons et calculer les ombrages
//Vec3 trace_ray(const Ray& ray, const std::vector<Sphere>& spheres, const std::vector<Light>& lights) {
//    auto intersection = intersect_scene(ray, spheres);
//    if (!intersection) {
//        return { 0.0f, 0.0f, 0.0f };  // Couleur de fond
//    }
//
//    Vec3 color = { 0.0f, 0.0f, 0.0f };
//    for (const auto& light : lights) {
//        Vec3 to_light = light.origin - intersection->point;
//        float light_distance = length(to_light);
//        Vec3 light_dir = normalize(to_light);
//
//        float cos_theta = std::max(0.0f, dot(intersection->normal, light_dir));
//        Vec3 intensity = intersection->albedo * (cos_theta / light_distance) * light.emission;
//
//        color = color + intensity;
//    }
//
//    return tonemap(color, 1.0f);
//}


//// Fonction principale pour lancer les rayons et calculer les ombrages
//Vec3 trace_ray(const Ray& ray, const Scene& scene, const std::vector<Light>& lights) {
//    auto intersection = intersect_scene(ray, scene);
//    if (!intersection) {
//        return { 0.0f, 0.0f, 0.0f };  // Couleur de fond
//    }
//
//    Vec3 color = { 0.0f, 0.0f, 0.0f };
//    for (const auto& light : lights) {
//        Vec3 to_light = light.origin - intersection->point;
//        float light_distance = length(to_light);
//        Vec3 light_dir = normalize(to_light);
//
//        float cos_theta = std::max(0.0f, dot(intersection->normal, light_dir));
//        Vec3 intensity = intersection->albedo * (cos_theta / light_distance) * light.emission;
//
//        color = color + intensity;
//    }
//
//    return tonemap(color, 1.0f);
//}




// Fonction principale pour lancer les rayons et calculer les ombrages
Vec3 trace_ray(const Ray& ray, const Scene& scene, const std::vector<Light>& lights) {
    auto intersection = scene.intersect(ray);
    if (!intersection) {
        return { 0.0f, 0.0f, 0.0f };
    }

    Vec3 color = { 0.0f, 0.0f, 0.0f };
    for (const auto& light : lights) {
        Vec3 to_light = light.origin - intersection->point;
        float light_distance = length(to_light);
        Vec3 light_dir = normalize(to_light);

        float cos_theta = std::max(0.0f, dot(intersection->normal, light_dir));
        Vec3 intensity = intersection->albedo * (cos_theta / light_distance) * light.emission;

        color = color + intensity;
    }

    return tonemap(color, 1.0f);
}


// Fonction pour enregistrer l'image
void save_image(const std::vector<uint8_t>& image, int width, int height, const std::string& filename) {
    std::ofstream file(filename, std::ios::out | std::ios::binary);
    file << "P6\n" << width << " " << height << "\n255\n";
    file.write(reinterpret_cast<const char*>(image.data()), image.size());
    file.close();
}

int main() {
    int width = 800;
    int height = 600;

    std::vector<Sphere> spheres = {
        {{0.0f, 0.0f, 200.0f}, 180.0f, {1.0f, 1.0f, 1.0f}},
        {{-300.0f, -300.0f, 200.0f}, 180.0f, {0.0f, 0.5f, 1.5f}},
        {{160.0f, 0.0f, 50.0f}, 50.0f, {2.0f, 2.0f, 2.0f}},
        {{0.0f, 50000.0f + 800.0f, 0.0f}, 50000.0f, {1.0f, 1.0f, 1.0f}}
    };

    std::vector<Light> lights = {
        {{5000.0f, 0.0f, 0.0f}, {400000.0f, 400000.0f, 400000.0f}},
        {{0.0f, -1000.0f, 0.0f}, {100000.0f, 0.0f, 0.0f}},
        {{-1000.0f, 1000.0f, 0.0f}, {0.0f, 100000.0f, 0.0f}}
    };

    Scene scene(spheres);



    float focal = 1000.0f;
    std::vector<uint8_t> image(width * height * 3);

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            Vec3 pixel = { float(x) * 2.0f - width, float(y) * 2.0f - height, 0.0f };
            Vec3 focal_point = { 0.0f, 0.0f, -focal };
            Vec3 direction = normalize(pixel - focal_point);

            Ray ray = { pixel, direction };
            Vec3 color = trace_ray(ray, scene,lights);

            image[(y * width + x) * 3 + 0] = static_cast<uint8_t>(color.x);
            image[(y * width + x) * 3 + 1] = static_cast<uint8_t>(color.y);
            image[(y * width + x) * 3 + 2] = static_cast<uint8_t>(color.z);
        }
    }

    save_image(image, width, height, "result.ppm");
    std::cout << "Image saved as result.ppm\n";
    return 0;
}

