#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "triangulation.h"
#include "texture.h"

using namespace std;

namespace CMU462 {

uint8_t* supersample_target;
Sampler2DImp sampler(SampleMethod method = TRILINEAR);

// Implements SoftwareRenderer //

void SoftwareRendererImp::draw_svg( SVG& svg ) {

  // allocate memory for supersample target
  supersample_target = new uint8_t[4 * sample_rate * target_w * sample_rate * target_h]();

  // draw all elements
  for ( size_t i = 0; i < svg.elements.size(); ++i ) {
    // set top level transformation
    transformation = canvas_to_screen;
    draw_element(svg.elements[i]);
  }

  // set top level transformation
  transformation = canvas_to_screen;

  // draw canvas outline
  Vector2D a = transform(Vector2D(    0    ,     0    )); a.x--; a.y++;
  Vector2D b = transform(Vector2D(svg.width,     0    )); b.x++; b.y++;
  Vector2D c = transform(Vector2D(    0    ,svg.height)); c.x--; c.y--;
  Vector2D d = transform(Vector2D(svg.width,svg.height)); d.x++; d.y--;

  // Disable outer polygon
  // rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
  // rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
  // rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
  // rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

  // resolve and send to render target
  resolve();
  
  // free the supersample target
  delete[] supersample_target;

}

void SoftwareRendererImp::set_sample_rate( size_t sample_rate ) {

  // Task 3: 
  // You may want to modify this for supersampling support
  this->sample_rate = sample_rate;

}

void SoftwareRendererImp::set_render_target( unsigned char* render_target,
                                             size_t width, size_t height ) {

  // Task 3: 
  // You may want to modify this for supersampling support
  this->render_target = render_target;
  this->target_w = width;
  this->target_h = height;

}

void SoftwareRendererImp::draw_element( SVGElement* element ) {

  // Task 4 (part 1):
  // Modify this to implement the transformation stack
  Matrix3x3 prev_transformation = transformation;
  transformation = transformation * element->transform;

  switch(element->type) {
    case POINT:
      draw_point(static_cast<Point&>(*element));
      break;
    case LINE:
      draw_line(static_cast<Line&>(*element));
      break;
    case POLYLINE:
      draw_polyline(static_cast<Polyline&>(*element));
      break;
    case RECT:
      draw_rect(static_cast<Rect&>(*element));
      break;
    case POLYGON:
      draw_polygon(static_cast<Polygon&>(*element));
      break;
    case FRACTAL:
      draw_fractal(static_cast<Fractal&>(*element));
      break;
    case ELLIPSE:
      draw_ellipse(static_cast<Ellipse&>(*element));
      break;
    case IMAGE:
      draw_image(static_cast<Image&>(*element));
      break;
    case GROUP:
      draw_group(static_cast<Group&>(*element));
      break;
    default:
      break;
  }

  transformation = prev_transformation;

}


// Primitive Drawing //

void SoftwareRendererImp::draw_point( Point& point ) {

  Vector2D p = transform(point.position);
  
  // fill the entire box to prevent shrinking during supersampling
  for (float x = floor(p.x) + (0.5f / sample_rate); x < floor(p.x) + 1.0f; x += (1.0f / sample_rate)) {
    for (float y = floor(p.y) + (0.5f / sample_rate); y < floor(p.y) + 1.0f; y += (1.0f / sample_rate)) {
      rasterize_point(x, y, point.style.fillColor);
    }
  }

}

void SoftwareRendererImp::draw_line( Line& line ) { 

  Vector2D p0 = transform(line.from);
  Vector2D p1 = transform(line.to);
  rasterize_line( p0.x, p0.y, p1.x, p1.y, line.style.strokeColor );

}

void SoftwareRendererImp::draw_polyline( Polyline& polyline ) {

  Color c = polyline.style.strokeColor;

  if( c.a != 0 ) {
    int nPoints = polyline.points.size();
    for( int i = 0; i < nPoints - 1; i++ ) {
      Vector2D p0 = transform(polyline.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polyline.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_rect( Rect& rect ) {

  Color c;
  
  // draw as two triangles
  float x = rect.position.x;
  float y = rect.position.y;
  float w = rect.dimension.x;
  float h = rect.dimension.y;

  Vector2D p0 = transform(Vector2D(   x   ,   y   ));
  Vector2D p1 = transform(Vector2D( x + w ,   y   ));
  Vector2D p2 = transform(Vector2D(   x   , y + h ));
  Vector2D p3 = transform(Vector2D( x + w , y + h ));
  
  // draw fill
  c = rect.style.fillColor;
  if (c.a != 0 ) {
    rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    rasterize_triangle( p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c );
  }

  // draw outline
  c = rect.style.strokeColor;
  if( c.a != 0 ) {
    rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    rasterize_line( p1.x, p1.y, p3.x, p3.y, c );
    rasterize_line( p3.x, p3.y, p2.x, p2.y, c );
    rasterize_line( p2.x, p2.y, p0.x, p0.y, c );
  }

}

void SoftwareRendererImp::draw_polygon( Polygon& polygon ) {

  Color c;

  // draw fill
  c = polygon.style.fillColor;
  if( c.a != 0 ) {

    // triangulate
    vector<Vector2D> triangles;
    triangulate( polygon, triangles );

    // draw as triangles
    for (size_t i = 0; i < triangles.size(); i += 3) {
      Vector2D p0 = transform(triangles[i + 0]);
      Vector2D p1 = transform(triangles[i + 1]);
      Vector2D p2 = transform(triangles[i + 2]);
      rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    }
  }

  // draw outline
  c = polygon.style.strokeColor;
  if( c.a != 0 ) {
    int nPoints = polygon.points.size();
    for( int i = 0; i < nPoints; i++ ) {
      Vector2D p0 = transform(polygon.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polygon.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}
  
void SoftwareRendererImp::draw_fractal( Fractal& fractal ) {
  // First print the seed point
  Vector2D point = fractal.seed;
  Vector2D transformed = transform(point);
  rasterize_point( transformed.x, transformed.y, fractal.style.fillColor );
  
  // Repeat for given number of iterations
  for ( long iter = 0; iter < fractal.iterations; iter++ ) {
    double rand = (double) std::rand() / RAND_MAX;
    double cumulative = 0;
    for ( std::vector<Transformation>::iterator tIt = fractal.transformations.begin();
            tIt != fractal.transformations.end(); tIt++) {
      cumulative += tIt->probability;
      
      if (rand <= cumulative ) {
        // Apply this selected transformation
        double new_x = tIt->a * point.x + tIt->b * point.y + tIt->e;
        double new_y = tIt->c * point.x + tIt->d * point.y + tIt->f;
        point = Vector2D( new_x, new_y );
        transformed = transform(point);
        rasterize_point( transformed.x, transformed.y, fractal.style.fillColor );
        break;
      }
  }
  }
}

void SoftwareRendererImp::draw_ellipse( Ellipse& ellipse ) {

  // Extra credit 
  Vector2D left_point = ellipse.center - ellipse.radius;
  Vector2D t_left_point = transform(left_point);
  Vector2D right_point = ellipse.center + ellipse.radius;
  Vector2D t_right_point = transform(right_point);

  Vector2D t_center = (t_left_point + t_right_point) / 2.0f;
  Vector2D t_radius = (t_right_point - t_left_point) / 2.0f;

  // Iterate samples in the bounding box to check coverage
  for ( float x = t_left_point.x + (0.5f / sample_rate); x < t_right_point.x + 1.0f; x += 1.0f / sample_rate ) {
    for ( float y = t_left_point.y + (0.5f / sample_rate); y < t_right_point.y + 1.0f; y += 1.0f / sample_rate ) {
      if ( ((x - t_center.x) * (x - t_center.x)) / (t_radius.x * t_radius.x)
           + ((y - t_center.y) * (y - t_center.y)) / (t_radius.y * t_radius.y) <= 1.0f ) {
        rasterize_point(x, y, ellipse.style.fillColor);
      }
    }
  }

}

void SoftwareRendererImp::draw_image( Image& image ) {

  Vector2D p0 = transform(image.position);
  Vector2D p1 = transform(image.position + image.dimension);

  rasterize_image( p0.x, p0.y, p1.x, p1.y, image.tex );
}

void SoftwareRendererImp::draw_group( Group& group ) {

  for ( size_t i = 0; i < group.elements.size(); ++i ) {
    draw_element(group.elements[i]);
  }

}

inline void uint8_to_float( float dst[4], unsigned char* src ) {
  uint8_t* src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
  dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8( unsigned char* dst, float src[4] ) {
  uint8_t* dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[0])));
  dst_uint8[1] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[1])));
  dst_uint8[2] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[2])));
  dst_uint8[3] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[3])));
}

// Rasterization //

// The input arguments in the rasterization functions 
// below are all defined in screen space coordinates

void SoftwareRendererImp::rasterize_point( float x, float y, Color element ) {

  // fill in the nearest pixel
  int sx = (int) floor(x * sample_rate);
  int sy = (int) floor(y * sample_rate);

  // check bounds
  if ( sx < 0 || sx >= sample_rate * target_w ) return;
  if ( sy < 0 || sy >= sample_rate * target_h ) return;

  // perform alpha blending
  Color canvas;
  uint8_to_float(&canvas.r, &supersample_target[4 * (sx + sy * sample_rate * target_w)]);
  Color blended;
  blended.a = 1.0f - (1.0f - element.a) * (1.0f - canvas.a);
  blended.r = (1.0f - element.a) * canvas.r + element.r;
  blended.g = (1.0f - element.a) * canvas.g + element.g;
  blended.b = (1.0f - element.a) * canvas.b + element.b;
  
  // fill sample
  float_to_uint8(&supersample_target[4 * (sx + sy * sample_rate * target_w)], &blended.r);

}

int get_octant(float x0, float y0, float x1, float y1) {
  float slope = (y1 - y0) / (x1 - x0);
  if (slope >= 0 && slope <= 1 && x0 < x1) {
    return 0;
  } else if (slope > 1 && y0 < y1) {
    return 1;
  } else if (slope < -1 && y0 < y1) {
    return 2;
  } else if (slope <= 0 && slope >= -1 && x1 < x0) {
    return 3;
  } else if (slope > 0 && slope <= 1 && x1 < x0) {
    return 4;
  } else if (slope > 1 && y1 < y0) {
    return 5;
  } else if (slope < -1 && y1 < y0) {
    return 6;
  } else /* if (slope < 0 && slope >= -1 && x0 < x1) */ {
    return 7;
  }
}

void switch_to_octant_zero_from( int octant, float &x, float &y ) {
  switch (octant) {
    case 0:
      // do nothing
      break;
    case 1:
      swap(x, y);
      break;
    case 2:
      x = -x;
      swap(x, y);
      break;
    case 3:
      x = -x;
      break;
    case 4:
      x = -x;
      y = -y;
      break;
    case 5:
      x = -x;
      y = -y;
      swap(x, y);
      break;
    case 6:
      y = -y;
      swap(x, y);
      break;
    default /* case 7 */:
      y = -y;
      break;
  }
}

void switch_from_octant_zero_to( int octant, float &x, float &y ) {
  switch (octant) {
    case 0:
      // do nothing
      break;
    case 1:
      swap(x, y);
      break;
    case 2:
      y = -y;
      swap(x, y);
      break;
    case 3:
      x = -x;
      break;
    case 4:
      x = -x;
      y = -y;
      break;
    case 5:
      x = -x;
      y = -y;
      swap(x, y);
      break;
    case 6:
      x = -x;
      swap(x, y);
      break;
    default /* case 7 */:
      y = -y;
      break;
  }
}

// Implementation of Bresenham's line algorithm
// based on pseudocode from Wikipedia
void SoftwareRendererImp::rasterize_line( float x0, float y0,
                                          float x1, float y1,
                                          Color color ) {
  // find out which octant the line belongs to
  int octant = get_octant(x0, y0, x1, y1);

  // switch to octant 0
  switch_to_octant_zero_from(octant, x0, y0);
  switch_to_octant_zero_from(octant, x1, y1);

  // start drawing the line
  float dx = x1 - x0;
  float dy = y1 - y0;
  float d = 2 * dy - dx;

  // draw the end point
  float draw_x = x0;
  float draw_y = y0;
  switch_from_octant_zero_to(octant, draw_x, draw_y);
  // fill the entire box to prevent shrinking during supersampling
  for (float sx = floor(draw_x) + (0.5f / sample_rate); sx < floor(draw_x) + 1.0f; sx += (1.0f / sample_rate)) {
    for (float sy = floor(draw_y) + (0.5f / sample_rate); sy < floor(draw_y) + 1.0f; sy += (1.0f / sample_rate)) {
      rasterize_point(sx, sy, color);
    }
  }

  float y = y0;
  for (float x = x0 + 1; x <= x1; x++) {
    if (d > 0) {
      y++;

      draw_x = x;
      draw_y = y;
      switch_from_octant_zero_to(octant, draw_x, draw_y);
      rasterize_point(draw_x, draw_y, color);      

      d += (2 * dy - 2 * dx);
    } else {
      draw_x = x;
      draw_y = y;
      switch_from_octant_zero_to(octant, draw_x, draw_y);
      for (float sx = floor(draw_x) + (0.5f / sample_rate); sx < floor(draw_x) + 1.0f; sx += (1.0f / sample_rate)) {
        for (float sy = floor(draw_y) + (0.5f / sample_rate); sy < floor(draw_y) + 1.0f; sy += (1.0f / sample_rate)) {
          rasterize_point(sx, sy, color);
        }
      }

      d += (2 * dy);
    }
  }
}

struct edge {
  float x_coefficient;
  float y_coefficient;
  float constant;
};

struct edge edge_equation( float x0, float y0,
                           float x1, float y1,
                           float x2, float y2 ) {
  struct edge e;
  e.x_coefficient = y1 - y0;
  e.y_coefficient = x0 - x1;
  e.constant = y0 * (x1 - x0) - x0 * (y1 - y0);

  // Invert sign to make sure third point is inside the line
  if ( x2 * e.x_coefficient + y2 * e.y_coefficient + e.constant > 0 ) {
    e.x_coefficient = -e.x_coefficient;
    e.y_coefficient = -e.y_coefficient;
    e.constant = -e.constant; 
  }

  return e;
}

void SoftwareRendererImp::rasterize_triangle( float x0, float y0,
                                              float x1, float y1,
                                              float x2, float y2,
                                              Color color ) {
  // Compute edge equations
  struct edge e[3];
  e[0] = edge_equation(x0, y0, x1, y1, x2, y2);
  e[1] = edge_equation(x1, y1, x2, y2, x0, y0);
  e[2] = edge_equation(x2, y2, x0, y0, x1, y1);

  // Compute top left and bottom right corners of the triangle's bounding box
  float top_left_x = floor(min(x0, min(x1, x2)));
  float top_left_y = floor(min(y0, min(y1, y2)));
  float bottom_right_x = floor(max(x0, max(x1, x2)));
  float bottom_right_y = floor(max(y0, max(y1, y2)));

  // Iterate samples in the bounding box to check coverage
  for ( float x = top_left_x + (0.5f / sample_rate); x < bottom_right_x + 1.0f; x += 1.0f / sample_rate ) {
    for ( float y = top_left_y + (0.5f / sample_rate); y < bottom_right_y + 1.0f; y += 1.0f / sample_rate ) {
      int n_inside_edges = 0;
      for ( int i = 0; i < 3; i++ ) {
        if ( x * e[i].x_coefficient + y * e[i].y_coefficient + e[i].constant <= 0 ) {
          n_inside_edges++;
        } else {
          break;
        }
      }

      if ( n_inside_edges == 3 ) {
        rasterize_point(x, y, color);
      }
    }
  }
}

void SoftwareRendererImp::rasterize_image( float x0, float y0,
                                           float x1, float y1,
                                           Texture& tex ) {
  // Task ?: 
  // Implement image rasterization
  // Compute transform from (x, y) space to (u, v) space
  float u_scale = tex.width / (x1 - x0);
  float v_scale = tex.height / (y1 - y0);

  // Fill texture
  float top_left_x = ceil(x0);
  float top_left_y = ceil(y0);
  float bottom_right_x = floor(x1);
  float bottom_right_y = floor(y1);

  // Iterate samples in the bounding box to check coverage
  for ( float x = top_left_x + (0.5f / sample_rate); x < bottom_right_x + 1.0f; x += 1.0f / sample_rate ) {
    for ( float y = top_left_y + (0.5f / sample_rate); y < bottom_right_y + 1.0f; y += 1.0f / sample_rate ) {
      if (x >= x0 && x <= x1 && y >= y0 && y<= y1) {
        float u = (x - x0) / (x1 - x0);
        float v = (y - y0) / (y1 - y0);
        Color c = sampler->sample_trilinear( tex, u, v, u_scale, v_scale );
        rasterize_point(x, y, c);
      }
    }
  }
}

// resolve samples to render target
void SoftwareRendererImp::resolve( void ) {

  // Task 3: 
  // Implement supersampling
  // You may also need to modify other functions marked with "Task 3".
  for ( int sx = 0; sx < target_w; sx++ ) {
    for ( int sy = 0; sy < target_h; sy++ ) {
      uint16_t box_sum[4] = { 0 };
      for ( int bx = sx * sample_rate; bx < (sx + 1) * sample_rate; bx++ ) {
        for ( int by = sy * sample_rate; by < (sy + 1) * sample_rate; by++ ) {
          box_sum[0] += supersample_target[4 * (bx + by * sample_rate * target_w)    ];
          box_sum[1] += supersample_target[4 * (bx + by * sample_rate * target_w) + 1];
          box_sum[2] += supersample_target[4 * (bx + by * sample_rate * target_w) + 2];
          box_sum[3] += supersample_target[4 * (bx + by * sample_rate * target_w) + 3];
        }
      }
      render_target[4 * (sx + sy * target_w)    ] = (uint8_t) (box_sum[0] / (sample_rate * sample_rate));
      render_target[4 * (sx + sy * target_w) + 1] = (uint8_t) (box_sum[1] / (sample_rate * sample_rate));
      render_target[4 * (sx + sy * target_w) + 2] = (uint8_t) (box_sum[2] / (sample_rate * sample_rate));
      render_target[4 * (sx + sy * target_w) + 3] = (uint8_t) (box_sum[3] / (sample_rate * sample_rate));
    }
  }

}


} // namespace CMU462
