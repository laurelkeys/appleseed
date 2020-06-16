
//
// This source file is part of appleseed.
// Visit https://appleseedhq.net/ for additional information and resources.
//
// This software is released under the MIT license.
//
// Copyright (c) 2020 Tiago Chaves, The appleseedhq Organization
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//

// Interface header.
#include "upsampleapplier.h"

// appleseed.renderer headers.
#include "renderer/utility/rgbcolorsampling.h"

// appleseed.foundation headers.
#include "foundation/image/canvasproperties.h"
#include "foundation/image/color.h"
#include "foundation/image/image.h"
#include "foundation/math/scalar.h"

using namespace foundation;

namespace renderer
{

//
// UpsampleApplier class implementation.
//

UpsampleApplier::UpsampleApplier(
    const UpsampleParams& params)
  : m_src_image(params.src_image)
  , m_src_width(params.src_image.properties().m_canvas_width)
  , m_src_height(params.src_image.properties().m_canvas_height)
{
}

void UpsampleApplier::release()
{
    delete this;
}

void UpsampleApplier::apply(
    Image&              image,
    const std::size_t   tile_x,
    const std::size_t   tile_y) const
{
    assert(tile_x < image.properties().m_tile_count_x);
    assert(tile_y < image.properties().m_tile_count_y);

    Tile& tile = image.tile(tile_x, tile_y);
    const std::size_t tile_width = tile.get_width();
    const std::size_t tile_height = tile.get_height();
    const Vector2u tile_offset(
        tile_x * image.properties().m_tile_width,
        tile_y * image.properties().m_tile_height);

    const std::size_t dst_width = image.properties().m_canvas_width;
    const std::size_t dst_height = image.properties().m_canvas_height;

    for (std::size_t y = 0; y < tile_height; ++y)
    {
        for (std::size_t x = 0; x < tile_width; ++x)
        {
            // Map the pixel coordinate from image to src_image.
            const float fx = static_cast<float>(x + tile_offset.x) / (dst_width - 1) * (m_src_width - 1);
            const float fy = static_cast<float>(y + tile_offset.y) / (dst_height - 1) * (m_src_height - 1);

            tile.set_pixel(x, y, dual_filter_upsample(m_src_image, fx, fy));
        }
    }
}

}   // namespace renderer
