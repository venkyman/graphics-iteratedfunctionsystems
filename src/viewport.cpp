#include "viewport.h"

#include "CMU462.h"

namespace CMU462 {

void ViewportImp::set_viewbox( float x, float y, float span ) {

  // Task 4 (part 2): 
  // Set svg to normalized device coordinate transformation. Your input
  // arguments are defined as SVG canvans coordinates.
  this->x = x;
  this->y = y;
  this->span = span;

  Matrix3x3 translate = Matrix3x3::identity();
  translate(0, 2) = - (x - span);
  translate(1, 2) = - (y - span);
  Matrix3x3 scale = Matrix3x3::identity();
  scale(0, 0) = 1.0f / (2.0f * span);
  scale(1, 1) = 1.0f / (2.0f * span);
  set_canvas_to_norm(scale * translate);
}

void ViewportImp::update_viewbox( float dx, float dy, float scale ) { 
  
  this->x -= dx;
  this->y -= dy;
  this->span *= scale;
  set_viewbox( x, y, span );
}

} // namespace CMU462
