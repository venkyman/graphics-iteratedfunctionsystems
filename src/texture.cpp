#include "texture.h"

#include <assert.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace CMU462 {

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

void Sampler2DImp::generate_mips(Texture& tex, int startLevel) {

  // NOTE(sky): 
  // The starter code allocates the mip levels and generates a level 
  // map simply fills each level with a color that differs from its
  // neighbours'. The reference solution uses trilinear filtering
  // and it will only work when you have mipmaps.

  // Task 7: Implement this

  // check start level
  if ( startLevel >= tex.mipmap.size() ) {
    std::cerr << "Invalid start level";
  }

  // allocate sublevels
  int baseWidth  = tex.mipmap[startLevel].width;
  int baseHeight = tex.mipmap[startLevel].height;
  int numSubLevels = (int)(log2f( (float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  tex.mipmap.resize(startLevel + numSubLevels + 1);

  int width  = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel& level = tex.mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width  = max( 1, width  / 2); assert(width  > 0);
    height = max( 1, height / 2); assert(height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(4 * width * height, 0);

  }

  // fill higher sublevels
  for(size_t i = 1; i < tex.mipmap.size(); ++i) {
    MipLevel& mip = tex.mipmap[i];
    MipLevel& lower_level = tex.mipmap[i - 1];

    for (size_t j = 0; j < 4 * mip.width * mip.height; j += 4) {
        int x = (j / 4) % mip.width;
        int y = (j / 4) / mip.width;

        Color top_left;
        uint8_to_float(&top_left.r, &lower_level.texels[4 * (2*x + 2*y * lower_level.width)]);
        Color top_right;
        uint8_to_float(&top_right.r, &lower_level.texels[4 * ((2*x + 1) + 2*y * lower_level.width)]);
        Color bottom_left;
        uint8_to_float(&bottom_left.r, &lower_level.texels[4 * (2*x + (2*y + 1) * lower_level.width)]);
        Color bottom_right;
        uint8_to_float(&bottom_right.r, &lower_level.texels[4 * ((2*x + 1) + (2*y + 1) * lower_level.width)]);
        
        Color avg = (0.25f * top_left + 0.25f * top_right + 0.25f * bottom_left + 0.25f * bottom_right);
        float_to_uint8( &mip.texels[j], &avg.r );
    }
  }

}

Color Sampler2DImp::sample_nearest(Texture& tex, 
                                   float u, float v, 
                                   int level) {

  // Task ?: Implement nearest neighbour interpolation
  if (level < tex.mipmap.size()) {
    MipLevel& miplevel = tex.mipmap[level];
    int iu = floor(u * miplevel.width);
    int iv = floor(v * miplevel.height);

    Color c;
    uint8_to_float(&c.r, &miplevel.texels[4 * (iu + iv * miplevel.width)]);
    return c;
  } else {
    // return magenta for invalid level
    return Color(1,0,1,1);
  }

}

Color Sampler2DImp::sample_bilinear(Texture& tex, 
                                    float u, float v, 
                                    int level) {
  
  // Task ?: Implement bilinear filtering
  if (level < tex.mipmap.size()) {
    MipLevel& miplevel = tex.mipmap[level];
    u = u * miplevel.width;
    v = v * miplevel.height;

    int lu;
    int uv;

    if (u >= floor(u) + 0.5f) {
      lu = floor(u);
    } else {
      lu = floor(u) - 1;
    }    

    if (v >= floor(v) + 0.5f) {
      uv = floor(v);
    } else {
      uv = floor(v) - 1;
    }

    float u_ratio = u - ((float) lu + 0.5f);
    float v_ratio = v - ((float) uv + 0.5f);

    Color top_left_color;
    if (lu >= 0 && uv >= 0) {
      uint8_to_float(&top_left_color.r, &miplevel.texels[4 * (lu + uv * miplevel.width)]);
    } else {
      top_left_color = Color(0, 0, 0, 1);
    }
    Color top_right_color;
    if (lu + 1 < miplevel.width && uv >= 0) {
      uint8_to_float(&top_right_color.r, &miplevel.texels[4 * (lu + 1 + uv * miplevel.width)]);
    } else {
      top_right_color = Color(0, 0, 0, 1);
    }
    Color bottom_left_color;
    if (lu >= 0 && uv + 1 < miplevel.height) {
      uint8_to_float(&bottom_left_color.r, &miplevel.texels[4 * (lu + (uv + 1) * miplevel.width)]);
    } else {
      bottom_left_color = Color(0, 0, 0, 1);
    }
    Color bottom_right_color;
    if (lu + 1 < miplevel.width && uv + 1 < miplevel.height) {
      uint8_to_float(&bottom_right_color.r, &miplevel.texels[4 * (lu + 1 + (uv + 1) * miplevel.width)]);
    } else {
      bottom_right_color = Color(0, 0, 0, 1);
    }

    return ((top_left_color * (1.0f - u_ratio) + top_right_color * u_ratio) * (1.0f - v_ratio)) +
            ((bottom_left_color * (1.0f - u_ratio) + bottom_right_color * u_ratio) * v_ratio);
  } else {
    // return magenta for invalid level
    return Color(1,0,1,1);
  }

}

Color Sampler2DImp::sample_trilinear(Texture& tex, 
                                     float u, float v, 
                                     float u_scale, float v_scale) {
  // Task 8: Implement trilinear filtering
  // compute d
  float L = max( u_scale, v_scale );
  float d = log2f(L);
  if (d <= 0) d = 0;

  // interpolate
  Color lower = sample_bilinear(tex, u, v, floor(d));
  Color upper = sample_bilinear(tex, u, v, ceil(d));
  float ratio = d - (float) floor(d);
  
  return lower * (1.0f - ratio) + upper * ratio;
}

} // namespace CMU462
