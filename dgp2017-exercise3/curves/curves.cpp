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
  // time step for smoothing
  double epsilon = 0.01;
  double tol = 0.005 ; // tolerance pour la difference de longueure acceptable
  double epsilon2 = 0.0001 ; // ratio de déplacement des points lors de la correction de la longueur

  int it = 0 ;
  Scalar mean_ptx = 0;
  Scalar mean_pty = 0;
  Scalar dist_init = 0 ;
  double delta_x,delta_y,expand_x,expand_y ;



  void laplacianSmoothing() {
    // Curve Smoothing - centroid (this function should do one iteration of smoothing)

      auto M = num_points ; // pour alléger l'écriture

//---------------- Calcul de la position du point moyen (centre du cercle estimé) et de la longueur total de la courbe
      it++ ;
      if (it == 1){
          for (int i = 0; i < M; ++i){
              mean_ptx = mean_ptx + points(0,i) ;
              mean_pty = mean_pty + points(1,i) ;
          }
          mean_ptx = mean_ptx / M ;
          mean_pty = mean_pty / M ;

        dist_init = calculate_length() ;
      }

//---------------- Déplacement des points celon l'argorithme mentionné dans l'ex.
      // Premier point fait séparément
      points(0,0) = (1- epsilon)*points(0,0) + epsilon*((points(0,M-1)+points(0,1))/2) ;
      points(1,0) = (1- epsilon)*points(1,0) + epsilon*((points(1,M-1)+points(1,1))/2) ;

      for (int i = 1; i < M-1; ++i){
          points(0,i) = (1- epsilon)*points(0,i) + epsilon*((points(0,i-1)+points(0,i+1))/2) ;
          points(1,i) = (1- epsilon)*points(1,i) + epsilon*((points(1,i-1)+points(1,i+1))/2) ;
      }

      // Dernier point fait séparément
      points(0,M-1) = (1- epsilon)*points(0,M-1) + epsilon*((points(0,M-2)+points(0,0))/2) ;
      points(1,M-1) = (1- epsilon)*points(1,M-1) + epsilon*((points(1,M-2)+points(1,0))/2) ;

//-----------------  correction de la longueur de la courbe en déplaçant les point (écartement radial d'un montant epsilon2 par rapport au centre pré_calculé)
      Scalar dist = 0 ;
      while (abs(dist - dist_init) > tol){
          dist = calculate_length() ;
cout << "dist:" << dist << endl ;

          if(dist<dist_init){
              for (int i = 0; i < M; ++i){
                  delta_x = points(0,i)- mean_ptx ;
                  delta_y = points(1,i)- mean_pty ;
                  expand_x = epsilon2/(sqrt( pow(delta_x,2) + pow(delta_y,2) )) * delta_x ;
                  expand_y = epsilon2/(sqrt( pow(delta_x,2) + pow(delta_y,2) )) * delta_y ;
                  points(0,i) = points(0,i) + expand_x ;
                  points(1,i) = points(1,i) + expand_y ;
              }
          }
          if(dist>dist_init){ // En principe on est pas plus petit mais c'est au cas ou la correction augmenterai trop la taille, on laisse ou pas ? avec les bon coef de correction je crois que c'est pas nécessaire
              for (int i = 0; i < M; ++i){
                  delta_x = points(0,i)- mean_ptx ;
                  delta_y = points(1,i)- mean_pty ;
                  expand_x = epsilon2/(sqrt( pow(delta_x,2) + pow(delta_y,2) )) * delta_x ;
                  expand_y = epsilon2/(sqrt( pow(delta_x,2) + pow(delta_y,2) )) * delta_y ;
                  points(0,i) = points(0,i) - expand_x ;
                  points(1,i) = points(1,i) - expand_y ;
              }
          }
      }
  }

  double calculate_length(){
      double length = 0 ;
      for (int i = 0; i < num_points-1; ++i){
          length = length + sqrt(pow(points(0,i)- points(0,i+1),2) + pow(points(1,i)-points(1,i+1),2)) ;
      }
      length = length + sqrt(pow(points(0,num_points-1)- points(0,0),2) + pow(points(1,num_points-1)-points(1,0),2)) ;
      return length;
  }

double mid_1_x,mid_1_y,mid_2_x,mid_2_y,perp_1_x,perp_1_y,perp_2_x,perp_2_y,a1,a2,b1,b2,norm_2 ;
double* center;
  void osculatingCircle() {
    // Curve Smoothing - osculating circle (again, this function should do one iteration of smoothing)
      auto M = num_points ; // pour alléger l'écriture
// Dans cette partie il reste a corriger les lignes ci-dessous car il y a une faute dans mon algorithme et je sais pas ou
// Ensuite il faut encore intésrer un calcul de la longueur de la courbe initial et un processus itératif de correction.
// Ceci peut être fait comme pour la partie du dessus a mon avis, on pourrais donc créer des fonction et les apeller dans les deux codes.
// j'espère que j'ai été clair ^^, si jamais redites moi je vous donnerai plus d'explication :D
      center = circumscribed_center(M-1,0,1) ;
      norm_2 = pow(center[0]-points(0,0),2) + pow(center[1]-points(1,0),2) ; // carré de la norm de la différence entre centre et le pt
      points(0,0) = points(0,0) + epsilon*(center[0]-points(0,0)/norm_2) ;
      points(1,0) = points(1,0) + epsilon*(center[1]-points(1,0)/norm_2)  ;

      for (int i = 1; i < M-1; ++i){
          center = circumscribed_center(i-1,i,i+1) ;
          norm_2 = pow(center[0]-points(0,i),2) + pow(center[1]-points(1,i),2) ; // carré de la norm de la différence entre centre et le pt
          points(0,i) = points(0,i) + epsilon*(center[0]-points(0,i)/norm_2) ;
          points(1,i) = points(1,i) + epsilon*(center[1]-points(1,i)/norm_2)  ;
//cout<<points(0,i)<<','<<points(1,i)<<endl ;
      }

      center = circumscribed_center(M-2,M-1,0) ;
      norm_2 = pow(center[0]-points(0,M-1),2) + pow(center[1]-points(1,M-1),2) ; // carré de la norm de la différence entre centre et le pt
      points(0,M-1) = points(0,M-1) + epsilon*(center[0]-points(0,M-1)/norm_2) ;
      points(1,M-1) = points(1,M-1) + epsilon*(center[1]-points(1,M-1)/norm_2)  ;

  }

  double* circumscribed_center(int i1,int i2,int i3){
      // Dans cette partie je calcul le point milieu des deux segments et le vecteur perpendiculaire a ces segments puis je crée les deux
      // equations des deux droites formées par ces vecteurs directeur ainsi que le point milieu et je cherche leur intersection afin de
      // trouver le centre du cercle !

      double centerr[2] ;
cout << points(0,i1) << ',' << points(0,i2) << ',' << points(0,i3) << endl ;
cout << points(1,i1) << ',' << points(1,i2) << ',' << points(1,i3) << endl ;
      mid_1_x = points(0,i1) + (points(0,i2) - points(0,i1))/2 ;
      mid_1_y = points(1,i1) + (points(1,i2) - points(1,i1))/2 ;

      mid_2_x = points(0,i2) + (points(0,i3) - points(0,i2))/2 ;
      mid_2_y = points(1,i2) + (points(1,i3) - points(1,i2))/2 ;

      // cas de la pente infinie ou nulle ? -> à voir comment on peut le résoudre
      perp_1_y = points(0,i2) - points(0,i1) ; // on prends l'inverse des coord pour obtenir le vecteur perpendiculaire
      perp_1_x = points(1,i2) - points(1,i1) ;

      perp_2_y = points(0,i3) - points(0,i2) ;
      perp_2_x = points(1,i3) - points(1,i2) ;

      a1 = (perp_1_y/perp_1_x) ;
      a2 = (perp_2_y/perp_2_x) ;

      b1 = mid_1_y - a1* mid_1_x ;
      b2 = mid_2_y - a2* mid_2_x ;

      centerr[0] = (b1-b2)/(a2-a1) ;
      centerr[1] = a1*centerr[0] + b1 ;
cout << centerr[0] << ',' << centerr[1] << endl ;
      return centerr ;
  }

// ============================================================================
// END OF Exercise 2 (do not thouch the rest of the code)
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
