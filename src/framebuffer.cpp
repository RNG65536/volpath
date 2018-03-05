#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <cstdint>
#include "vector.hpp"
#include "framebuffer.hpp"
#include "utils.hpp"

using std::cout;
using std::endl;

static void dumpHDR(const char *filename, const float *m_buffer, int m_width, int m_height)
{
    struct HDRPixel
    {
        uint8_t r, g, b, e;
        HDRPixel() = default;
        HDRPixel(uint8_t r, uint8_t g, uint8_t b, uint8_t e) : r(r), g(g), b(b), e(e)
        {

        }
        HDRPixel(const vec3& c)
        {
            float d = f_max(f_max(c.x, c.y), c.z);
            if (d <= 1e-32f)
            {
                r = g = b = e = 0;
                return;
            }
            int _e;
            float m = frexp(d, &_e); // d = m * 2^e
            d = m * 256.0f / d;
            r = uint8_t(c.x * d);
            g = uint8_t(c.y * d);
            b = uint8_t(c.z * d);
            e = uint8_t(_e + 128);
        }

        uint8_t operator[](int i) const
        {
            return (&r)[i];
        }
    };

    std::ofstream ofs(filename, std::ios::out | std::ios::binary);
    if (!ofs.is_open())
    {
        cout << "cannot write to " << filename << endl;
    }

    ofs << "#?RADIANCE" << endl;
    ofs << "# Made with custom writer" << endl;
    ofs << "FORMAT=32-bit_rle_rgbe" << endl;
    ofs << "EXPOSURE=1.0" << endl;
    ofs << endl;

    ofs << "-Y " << m_height << " +X " << m_width << endl;

    for (int j = m_height - 1; j >= 0; --j)
    {
        std::vector<HDRPixel> line(m_width);
        for (int i = 0; i < m_width; i++)
        {
            int offset = (i + j * m_width) * 3;
            HDRPixel p = HDRPixel(vec3(m_buffer[offset], m_buffer[offset + 1], m_buffer[offset + 2]));
            line[i] = p;
        }
        ofs << uint8_t(2) << uint8_t(2);
        ofs << uint8_t((m_width >> 8) & 0xFF)
            << uint8_t(m_width & 0xFF);
        for (int k = 0; k < 4; k++)
        {
            for (int cursor = 0; cursor < m_width;)
            {
                const int cursor_move = std::min(127, m_width - cursor);
                ofs << uint8_t(cursor_move);
                for (int i = cursor; i < cursor + cursor_move; i++)
                {
                    ofs << uint8_t(line[i][k]);
                }
                cursor += cursor_move;
            }
        }
    }
    ofs.close();
}


static void save_hdr_image(const float *fb, int width, int height)
{
    dumpHDR("_additional.hdr", fb, width, height);

    struct ppm_writer
    {
        ppm_writer(){}
        float clamp(float x){ return x < 0 ? 0 : x>1 ? 1 : x; }
        void operator()(const float *fb, int width, int height, const char *fn)
        {
            std::ofstream ofs(fn, std::ios::out | std::ios::binary);
            ofs << "P6\n" << width << " " << height << "\n255\n";
            unsigned char *bufByte = new unsigned char[width * height * 3];
            int n = 0;
            for (int j = height - 1; j >= 0; j--)
            {
                for (int i = 0; i < width; i++)
                {
                    int offset = (i + j * width) * 3;
//                     bufByte[n++] = (unsigned char)(std::pow(clamp(fb[offset + 0]), 1.0f / 2.2f) * 255);
//                     bufByte[n++] = (unsigned char)(std::pow(clamp(fb[offset + 1]), 1.0f / 2.2f) * 255);
//                     bufByte[n++] = (unsigned char)(std::pow(clamp(fb[offset + 2]), 1.0f / 2.2f) * 255);
                    bufByte[n++] = (unsigned char)(clamp(fb[offset + 0]) * 255);
                    bufByte[n++] = (unsigned char)(clamp(fb[offset + 1]) * 255);
                    bufByte[n++] = (unsigned char)(clamp(fb[offset + 2]) * 255);
                }
            }
            ofs.write(reinterpret_cast<char*>(bufByte), width * height * 3);
            ofs.close();
            delete[] bufByte;
        }
    }
    save_ppm;

    save_ppm(&fb[0], width, height, "_additional.ppm");
}

void FrameBuffer::dump() const
{
    const float *fb1 = reinterpret_cast<const float*>(m_buffer.data());
    printf("color range: [%f, %f] \n",
        *std::min_element(fb1, fb1 + m_buffer.size() * 3),
        *std::max_element(fb1, fb1 + m_buffer.size() * 3));

    //not tonemapped
    save_hdr_image(fb1, m_width, m_height);
}

std::vector<vec3>& FrameBuffer::buffer()
{
    return m_buffer;
}

void FrameBuffer::scaleBrightness(float s)
{
    struct vec_scale_op{
        float _s;
        vec_scale_op(float s) :_s(s){}
        vec3 operator()(vec3 x){
            return x*_s;
        }
    };
    std::transform(m_buffer.begin(), m_buffer.end(), m_buffer.begin(), vec_scale_op(s));
}

FrameBuffer::FrameBuffer(int width, int height) : m_width(width), m_height(height)
{
    m_buffer.resize(m_width * m_height);
    std::fill(m_buffer.begin(), m_buffer.end(), vec3(0));
}
