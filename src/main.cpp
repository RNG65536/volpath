#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <cassert>
#include <string>
#include <cstdint>
using std::cout;
using std::endl;

#include "constants.h"
#include "vector.h"
#include "utils.h"
#include "framebuffer.h"
#include "light.h"
#include "ray.h"
#include "texture3d.h"
#include "scene.h"
#include "camera.h"

namespace P
{
    int width, height;
    int spp;
    int trace_depth;

    float fov;
    vec3 camera_origin;
    vec3 camera_direction;

    float density_scale;
    float HG_mean_cosine;

    bool show_bg;

    void param();
};

void P::param()
{
    float dist = 3.0f;
    show_bg = true;

    float rotate = 0.5f;
    camera_origin = vec3(dist*sinf(rotate), 1.2f, dist*cosf(rotate));
    camera_direction = vec3(0.0f);
    fov = 30.0f;

    width = 600;
    height = 600;
    spp = 5;
    trace_depth = 10;// 100;

    density_scale = 10.0f; //  50; //100; //400;
    HG_mean_cosine = 0.0f;// 0.7;//higher improves performance
}

class HeterogeneousMedium
{
    Texture3D *m_density = nullptr;
    vec3 m_sigma_s, m_sigma_a, m_sigma_t;
    vec3 m_albedo;
    vec3 m_sigma_t_max; // sigma is directly modulated by density
    vec3 m_inv_sigma_t_max; // sigma is directly modulated by density
    float m_density_max;
    float m_inv_density_max;

public:
    HeterogeneousMedium(Texture3D *vol, float scale)
    {
        m_density = vol;
        float max_density = vol->maxValue();
        printf("maximum density is %f\n", max_density);
        m_density_max = max_density;
        m_inv_density_max = 1.0f / m_density_max;

        m_sigma_s = vec3(0.70f, 1.22f, 1.90f) * scale;
        m_sigma_a = vec3(0.0014f, 0.0025f, 0.0142f) * scale;

        m_sigma_t = m_sigma_s + m_sigma_a;
        m_albedo = m_sigma_s / m_sigma_t;
        m_sigma_t_max = m_sigma_t * max_density;
        m_inv_sigma_t_max = vec3(1.0f) / m_sigma_t_max;
    }
    float sigmaT(const vec3& pos, int channel) const
    {
        return m_sigma_t[channel] * m_density->fetch(pos);
    }
    float albedo(int channel)
    {
        return m_albedo[channel];
    }
    Texture3D *density()
    {
        return m_density;
    }

    // continue only if null collision
    bool isNotNullCollision(const vec3& pos, float rnd) const
    {
//         return rnd < m_sigma_t.x * m_density->fetch(pos) * m_inv_sigma_t_max.x;
        return rnd < m_density->fetch(pos) * m_inv_density_max;
    }
    float sampleFreeDistance(float rnd, int channel) const
    {
        return -std::log(1.0f - rnd) * m_inv_sigma_t_max[channel];
    }
};

Light light;
Camera *cam = nullptr;
HeterogeneousMedium *medium = nullptr;

vec3 sampleHG(float g, float e1, float e2)
{
    float s = 2.0f * e1 - 1.0f, f = (1.0f - g * g) / (1.0f + g * s);
    float cost = 0.5f * (1.0f / g) * (1.0f + g * g - f * f), sint = std::sqrt(1.0f - cost * cost);
    return vec3(std::cos(2.0f * M_PI * e2) * sint, std::sin(2.0f * M_PI * e2) * sint, cost);
}

float HG_phase_function(float g, float cos_t){
    return (1.0f - g * g) / (4.0f * M_PI * std::pow(1.0f + g * g - 2 * g * cos_t, 1.5f));
}

float stochasticTransmittance(
    const vec3& start_point,
    const vec3& end_point,
    int channel)
{
    Ray light_ray(start_point, normalize(end_point - start_point));

    float t_near, t_far;
    bool shade_vol = intersect_vol(light_ray, medium->density()->min, medium->density()->max, t_near, t_far);
    if (!shade_vol)
    {
        return 1;
    }

    float max_t = f_min(t_far, distance(start_point, end_point));

    int nsamp = 2;
    float count = 0; // the samples that reached the destination by delta tracking

    for (int n = 0; n < nsamp; n++)
    {
        /// woodcock tracking
        float dist = t_near;

        while (true)
        {
            dist += medium->sampleFreeDistance(rand01(), channel);
            if (dist >= max_t)
            {
                count += 1;
                break;
            }
            vec3 pos = light_ray.at(dist);

            if (medium->isNotNullCollision(pos, rand01()))
            {
                break;
            }
        }
    }
    return count / nsamp;
}

float volumetricPathTacing(
    Ray& ray,
    const Texture3D& pdensity_volume,
    int channel)
{
    Ray cr(ray);
    float radiance = (0.0f);
    float throughput = (1.0f);

    /// for woodcock tracking
    int max_depth(P::trace_depth);
    for (int depth = 0; depth < max_depth; depth++)
    {
        float t_near, t_far;
        bool shade_vol = intersect_vol(cr, pdensity_volume.min, pdensity_volume.max, t_near, t_far);

        if (!shade_vol && depth == 0 && P::show_bg)
        {
            radiance = radiance + (light.Li(cr.o, cr.d)[channel] * throughput);
            break;
        }

        vec3 front = cr.at(t_near);
        vec3 back = cr.at(t_far);

        if (shade_vol)
        {
            /// woodcock tracking
            vec3 pos = front;//current incident radiance evaluation point
            float dist = t_near;

            int through = 0;
            while (true)
            {
                dist += medium->sampleFreeDistance(rand01(), channel);
                if (dist >= t_far){
                    through = 1;//transmitted through the volume, probability is 1-exp(-optical_thickness)
                    break;
                }
                pos = cr.at(dist);
                if (medium->isNotNullCollision(pos, rand01()))
                {
                    break;
                }
            }

            if (0 == through)
            {
                throughput = throughput * medium->albedo(channel); //all subsequent light evaluations are scattered

                // the sun light may be strongly directional, but is always scattered by atmosphere,
                // hence modeled by IBL (depends on Mie phase function for directionality)
                {
                    vec3 Li;
                    float pdfW;
                    vec3 lpos, ldir;
                    light.sample(pos, pdfW, lpos, ldir, Li);

                    float transmittance = stochasticTransmittance(pos, lpos, channel);
                    float attenuated_radiance = Li[channel] * transmittance;
                    radiance = radiance +
                        (attenuated_radiance * throughput) *
                        (HG_phase_function(P::HG_mean_cosine, dot(cr.d, ldir))
                        / pdfW);
                }

                vec3 dir = sampleHG(P::HG_mean_cosine, rand01(), rand01());
                vec3 ref_dir(2, 3, 5);
                ref_dir = normalize(ref_dir);
                vec3 u = normalize(cross(cr.d, ref_dir));
                vec3 v = cross(cr.d, u);

                dir = u * dir.x + v * dir.y + cr.d * dir.z;//by construction of the sample coordinates, dir is guaranteed to be unit
                cr = Ray(pos, dir);
            }
            else
            {
                if (depth == 0 && P::show_bg)
                {
                    radiance = radiance + (light.Li(cr.o, cr.d)[channel] * throughput);
                }

                break;
            }

        }
    }

    if (radiance != radiance)
    {
        radiance = 0;
    }

    return radiance;
}

void render(const Texture3D& vol)
{
    FrameBuffer framebuffer(P::width, P::height);

    for (int k = 0; k < P::spp; k++)
    {
        for (int j = 0; j < P::height; j++)
        {
            for (int i = 0; i < P::width; i++)
            {
                framebuffer.buffer()[i + j * P::width] += vec3(
                    volumetricPathTacing(cam->pixelRay(i, j), vol, 0),
                    volumetricPathTacing(cam->pixelRay(i, j), vol, 1),
                    volumetricPathTacing(cam->pixelRay(i, j), vol, 2)
                    );
            }
        }
        printf("\rfinished %.2f%%", 100.0f * float(k + 1) / float(P::spp));
    }

    framebuffer.scaleBrightness(1.0f / P::spp);
    framebuffer.dump();

    printf("\rrendering finished\n");
}

int main(int argc, char *argv[])
{
    P::param();
    cam = new Camera(P::camera_origin, P::camera_direction, P::fov, P::width, P::height);

    printf("DEBUG:: begin program\n");
    Texture3D vol(128, 1, vec3(-0.5));

    vol.initChecker();
//     vol.load_binary("smoke.vol");

    medium = new HeterogeneousMedium(&vol, P::density_scale);

    printf("init complete\n");
    printf("rendering\n");
    render(vol);

    return 0;
}