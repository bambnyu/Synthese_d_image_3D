#include <iostream>
#include <fstream>
#include <cmath>
#include <optional>
#include <vector>



// structure de Vecteur 3D
struct Vec3 {
    float x, y, z;

    // Operatuer de soustraction de vecteur
    Vec3 operator-(const Vec3& other) const {
        return { x - other.x, y - other.y, z - other.z };
    }
};

// structure de Sphère
struct Sphere {
    float radius;
    Vec3 center;
};

// structure de Rayon
struct Rayon {
    Vec3 origin;
    Vec3 direction;
};






// Produit scalaire
float dot(const Vec3& a, const Vec3& b) {
    return ((a.x * b.x) + (a.y * b.y) + (a.z * b.z));
}

// Calcul de la norme d'un vecteur 3D pour avoir sa longueur
float length(const Vec3& a) {
    return std::sqrt(dot(a, a));
}

// Norme au carré
float length_squarred(const Vec3& a) {
    return dot(a, a);
}

// calcul de la distance entre l'origine du rayon et le point d'intersection
float get_intersection_distance(const Rayon& ray, float t) {
    return t * length(ray.direction);
}




// calcul de l'intersection entre un rayon et une sphère
std::optional<float>  intersect_sphere(const Rayon& ray, const Sphere& sphere) {
    Vec3 oc = ray.origin - sphere.center; // vecteur entre l'origine du rayon et le centre de la sphère

    float a = length_squarred(ray.direction); // norme du vecteur direction du rayon
    float b = 2.0f * dot(oc, ray.direction); // produit scalaire entre le vecteur oc et le vecteur direction du rayon
    float c = length_squarred(oc) - sphere.radius * sphere.radius; // norme du vecteur oc - rayon de la sphère

    float delta = b * b - 4.0f * a * c;

    if (delta >= 0.0f) {
        // calcul des deux solutions de l'équation du second degré
        float t1 = (-b - std::sqrt(delta)) / (2.0f * a);
        float t2 = (-b + std::sqrt(delta)) / (2.0f * a);

        if (t1 >= 0.0f) {
            return t1;
        }
        else if (t2 >= 0.0f) {
            return t2;
        }
    }
   return std::nullopt;
}



// ecriture de l'image en format PPM
void save_image(const std::vector<unsigned char>& image, int width, int height, const std::string& filename) {
    std::ofstream file(filename, std::ios::out | std::ios::binary); // ouverture du fichier en mode binaire
    file << "P6\n" << width << " " << height << "\n255\n"; // entête du fichier PPM
    file.write(reinterpret_cast<const char*>(image.data()), image.size()); // écriture des données de l'image
    file.close();
}



int main() {
    int w = 800;
    int h = 600;

    std::vector<unsigned char> img(w * h * 3, 0); // initialisation de l'image *3 pour les 3 couleurs

    float radius = 180.0f;
    Sphere sphere = { radius, {0.0f, 0.0f, 200.0f} };
    float focal = 10000.0f;

    // la faire la boucle de remplissage et de test pour l'intersection cercle
    for (int py = 0; py < h; ++py) {
        float y = static_cast<float>(py);

        for (int px = 0; px < w; ++px) {
            float x = static_cast<float>(px);

            // Pixel as a 3D point
            Vec3 pixel = { x * 2.0f - w, y * 2.0f - h, 0.0f };
            Vec3 focal_point = { 0.0f, 0.0f, -focal };
            Vec3 direction = pixel - focal_point;

            Rayon ray = { pixel, direction };

            auto it = intersect_sphere(ray, sphere);

            if (it.has_value()) {
                float distance = get_intersection_distance(ray, it.value());
                unsigned char v = static_cast<unsigned char>(std::min(distance * 1.0f, 255.0f));

                // Set pixel color (grayscale)
                img[(py * w + px) * 3 + 0] = v; // Red
                img[(py * w + px) * 3 + 1] = v; // Green
                img[(py * w + px) * 3 + 2] = v; // Blue
            }
            else {
                // Set background color (e.g., dark blue-ish)
                img[(py * w + px) * 3 + 0] = 72;
                img[(py * w + px) * 3 + 1] = 61;
                img[(py * w + px) * 3 + 2] = 139;
            }
        }
    }


    save_image(img, w, h, "result.ppm"); // j'enregistre l'image

    std::cout << "Image enregistré comme result.ppm" << std::endl;
    return 0;
}

