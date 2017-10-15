#include <OpenGP/GL/TrackballWindow.h>
#include <OpenGP/SurfaceMesh/GL/SurfaceMeshRenderShaded.h>
#include <OpenGP/SurfaceMesh/GL/SurfaceMeshRenderFlat.h>
#include <OpenGP/GL/PointsRenderer.h>
#include <OpenGP/GL/SegmentsRenderer.h>
#include <random>

using namespace OpenGP;

using std::cout ;
using std::endl ;
struct MainWindow : public TrackballWindow {
    PointsRenderer render_points = PointsRenderer ();
    SegmentsRenderer render_segments = SegmentsRenderer ();

    MatMxN points;
    MatMxN points_3d_render;
    int num_points;
    double radius = 0.3;
    SegmentsRenderer::Segments segments;

    // ============================================================================
    // Exercise 2 : fill the 2 functions below (see PDF for instructions)
    // To test your implementation, use the S key for laplacian smoothing and the
    // C key for the osculating circle.
    // Hint : try to play with epsilon
    // ============================================================================
    Scalar epsilon = 0.001; // time step for smoothing
    Scalar tol = 0.005; // tolerance pour la difference de longueure acceptable
    Scalar epsilon2 = 0.0001; // ratio de déplacement des points lors de la correction de la longueur

    int it = 0 ;
    Scalar initial_length = 0.;

    void laplacianSmoothing() {
        // Curve Smoothing - centroid (this function should do one iteration of smoothing)
        auto m = num_points; // pour alléger l'écriture

        if (it++ == 0){
            initial_length = calculate_length();
        }

        // Déplacement des points celon l'algorithme mentionné dans l'ex.
        for (int i = 0; i < m - 1; ++i){
            points(0, i) = (1 - epsilon) * points(0, i) + epsilon * ((points(0, (m + i - 1) % m) + points(0, i + 1)) / 2.) ;
            points(1, i) = (1 - epsilon) * points(1, i) + epsilon * ((points(1, (m + i - 1) % m) + points(1, i + 1)) / 2.) ;
        }

        // Correction de la longueur de la courbe en déplaçant les point (écartement radial d'un montant epsilon2 par rapport au centre pré_calculé)
        correctLength();
    }


    void osculatingCircle() {
        // Curve Smoothing - osculating circle (again, this function should do one iteration of smoothing)
        auto m = num_points; // pour alléger l'écriture

        // Dans cette partie il reste a corriger les lignes ci-dessous car il y a une faute dans mon algorithme et je sais pas où
        // Ceci peut être fait comme pour la partie du dessus a mon avis, on pourrais donc créer des fonction et les apeller dans les deux codes.
        // j'espère que j'ai été clair ^^, si jamais redites moi je vous donnerai plus d'explication :D

        if (it++ == 0){
            initial_length = calculate_length();
        }

        Scalar diff_x, diff_y, norm_2 = 1.;
        Vec2 center;
        for (int i = 0; i < m - 1; i++) {
            center = circumscribedCenter((m + i - 1) % m, i, i + 1);
            diff_x = center[0] - points(0, i);
            diff_y = center[1] - points(1, i);
            norm_2 = diff_x * diff_x + diff_y * diff_y; // carré de la norme de la différence entre centre et le pt
            points(0, i) += epsilon * (diff_x / norm_2);
            points(1, i) += epsilon * (diff_y / norm_2);
        }

        correctLength();
    }

    Vec2 circumscribedCenter(int i1, int i2, int i3){
        // Dans cette partie je calcul le point milieu des deux segments et le vecteur perpendiculaire a ces segments puis je crée les deux
        // equations des deux droites formées par ces vecteurs directeur ainsi que le point milieu et je cherche leur intersection afin de
        // trouver le centre du cercle !
        Vec2 center;

        Vec2 perp_1 = {points(1,i2) - points(1,i1),
                       points(0,i1) - points(0,i2)};
        Vec2 perp_2 = {points(1,i3) - points(1,i2),
                       points(0,i2) - points(0,i3)};

        double a1, a2, b1, b2;

        // Géometriquement faux, mais permet de supporter les cas où les trois points seraient alignés
        // aka perp_1 et perp_2 sont verticales
        if (perp_1[0] == 0 && perp_2[0] == 0){
            center[0] = points(0, i3) - points(0, i1);
            center[1] = points(1, i2);
            return center;
        }

        Vec2 mid_1 = {points(0, i1) + (points(0, i2) - points(0, i1)) / 2.,
                      points(1, i1) + (points(1, i2) - points(1, i1)) / 2.};
        Vec2 mid_2 = {points(0, i2) + (points(0, i3) - points(0, i2)) / 2.,
                      points(1, i2) + (points(1, i3) - points(1, i2)) / 2.};

        // perp_1 est verticale
        if (perp_1[0] == 0) {
            a2 = perp_2[1] / perp_2[0];
            b2 = mid_2[1] - a2 * mid_2[0];
            center[0] = mid_1[0];
            center[1] = a2 * center[0] + b2;
            return center;
        }

        // perp_2 est verticale
        if (perp_2[0] == 0) {
            a1 = perp_1[1] / perp_1[0];
            b1 = mid_1[1] - a1 * mid_1[0];
            center[0] = mid_2[0];
            center[1] = a1 * center[0] + b1;
            return center;
        }

        // Cas classique, deux pentes non verticales
        a1 = perp_1[1] / perp_1[0];
        a2 = perp_2[1] / perp_2[0];
        b1 = mid_1[1] - a1 * mid_1[0];
        b2 = mid_2[1] - a2 * mid_2[0];

        center[0] = (b2 - b1) / (a1 - a2);
        center[1] = a1 * center[0] + b1;
        return center;

    }

    void correctLength(){
        Scalar length = 0.;
        Scalar delta_x, delta_y;
        Scalar expand_x, expand_y;
        Vec2 center = computeCenter();
        Scalar center_x = center[0], center_y = center[1];
        while (std::abs(length - initial_length) > tol){
            for (int i = 0; i < num_points - 1; ++i){
                delta_x = points(0,i) - center_x;
                delta_y = points(1,i) - center_y;
                expand_x = epsilon2 / (sqrt(delta_x * delta_x + delta_y * delta_y)) * delta_x;
                expand_y = epsilon2 / (sqrt(delta_x * delta_x + delta_y * delta_y)) * delta_y;

                if (length > initial_length) {
                    points(0,i) -= expand_x;
                    points(1,i) -= expand_y;
                } else {
                    points(0,i) += expand_x;
                    points(1,i) += expand_y;
                }
            }

            length = calculate_length();
        }
        cout << "initial length: " << initial_length << endl;
        cout << "length: " << length << endl;
    }

    Scalar calculate_length(){
        Scalar dx, dy, length = 0 ;
        for (int i = 0; i < num_points - 1; ++i){
            dx = points(0,i) - points(0, i + 1);
            dy = points(1,i) - points(1, i + 1);
            length += sqrt(dx * dx + dy * dy) ;
        }
        return length;
    }

    // Calcul de la position du point moyen (centre du cercle estimé) et de la longueur total de la courbe
    Vec2 computeCenter(){
        Scalar center_x = 0., center_y = 0.;
        for (int i = 0; i < num_points - 1; i++){
            center_x += points(0,i);
            center_y = points(1,i);
        }
        center_x /= num_points - 1;
        center_y /= num_points - 1;
        return {center_x, center_y};
    }

    // ============================================================================
    // END OF Exercise 2 (do not touch the rest of the code)
    // ============================================================================

    void generateRandomizedClosedPolyline() {
        std::default_random_engine generator;
        std::uniform_real_distribution<double> distribution(0., 5*3e-2);

        Vec2 center(3e-2, 2e-3);
        it = 0 ;

        points = MatMxN::Zero(2, num_points);
        for (int i = 0; i < num_points; ++i)
        {
            double frac = static_cast<double>(i) / static_cast<double>(num_points);
            points(0, i) = center(0) + radius * cos (2. * M_PI * frac) + distribution(generator);
            points(1, i) = center(1) + radius * sin (2. * M_PI * frac) + distribution(generator);
        }
    }

    void render () {

        // Prepare the render points
        points_3d_render = MatMxN::Zero(3, points.cols());
        points_3d_render.block(0, 0, 2, points.cols()) = points;

        // Rebuild the segments
        segments.clear();
        for (int i = 0; i < points_3d_render.cols(); ++i) {
            segments.push_back({ points_3d_render.col(i), points_3d_render.col((i+1) % points_3d_render.cols()) });
        }
        render_points.init_data(points_3d_render);
        render_segments.init_data(segments);
    }

    MainWindow(int argc, char** argv) : TrackballWindow("2D Viewer", 640, 480) {
        num_points = 30;
        generateRandomizedClosedPolyline();

        this->scene.add(render_points);
        this->scene.add(render_segments);

        render();
    }

    bool key_callback(int key, int scancode, int action, int mods) override {
        TrackballWindow::key_callback(key, scancode, action, mods);
        if (key == GLFW_KEY_S && action == GLFW_RELEASE)
        {
            laplacianSmoothing();
        }
        else if (key == GLFW_KEY_C && action == GLFW_RELEASE)
        {
            osculatingCircle();
        }
        else if (key == GLFW_KEY_1 && action == GLFW_RELEASE)
        {
            num_points = 30;
            radius = 0.3;
            generateRandomizedClosedPolyline();
        }
        else if (key == GLFW_KEY_2 && action == GLFW_RELEASE)
        {
            num_points = 50;
            radius = 0.1;
            generateRandomizedClosedPolyline();
        }
        else if (key == GLFW_KEY_3 && action == GLFW_RELEASE)
        {
            num_points = 100;
            radius = 0.1;
            generateRandomizedClosedPolyline();
        }
        else if (key == GLFW_KEY_4 && action == GLFW_RELEASE)
        {
            num_points = 150;
            radius = 0.1;
            generateRandomizedClosedPolyline();
        }

        render();
        return true;
    }
};


int main(int argc, char** argv)
{
    MainWindow window(argc, argv);
    return window.run();
}
