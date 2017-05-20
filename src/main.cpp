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
#include <memory>
#include <utility>

using std::cout;
using std::endl;

#include "constants.h"
#include "vector.h"
#include "utils.h"
#include "framebuffer.h"
#include "skydome.h"
#include "light.h"
#include "ray.h"
#include "texture3d.h"
#include "scene.h"
#include "camera.h"
#include "fractal.h"
#include "medium.h"

namespace P
{
    int width, height;
    int spp;
    int trace_depth;

    float fov;
    vec3 camera_origin;
    vec3 camera_lookat;

    float density_scale;
    float HG_mean_cosine;
    vec3 sigma_s, sigma_a;

    bool show_bg;

    void param();
};

void P::param()
{
    float dist = 3.0f;
    show_bg = true;

    float rotate = 0.5f;
//     camera_origin = vec3(dist*sinf(rotate), 1.2f, dist*cosf(rotate));
    camera_origin = vec3(dist*sinf(rotate), -0.55f, dist*cosf(rotate));
    camera_lookat = vec3(0.0f, 0.2f, 0.0f);
    fov = 30.0f;

    width = 600;
    height = 600;
    spp = 100;
    trace_depth = 500;

    // sigma_t_prime = sigma_a + sigma_s * (1 - g)
    // for scattering dominated media, should scale with 1 / (1 - g) to approximate appearance
    density_scale = 1000.0f; //  20.0f;
    HG_mean_cosine = 0.877f;// 0.7f;// 

//     sigma_s = vec3(0.70f, 1.22f, 1.90f) * density_scale;
//     sigma_a = vec3(0.0014f, 0.0025f, 0.0142f) * density_scale;
    sigma_s = vec3(1.0f, 1.0f, 1.0f) * density_scale;
    sigma_a = vec3(0.0001f, 0.0001f, 0.0001f) * density_scale;
}

Light light;
std::unique_ptr<Camera> cam;
std::unique_ptr<HGPhaseFunction> phase;
std::unique_ptr<HeterogeneousMedium> medium;
std::unique_ptr<Texture3D> vol;

float stochasticTransmittance(
    const vec3& start_point,
    const vec3& end_point,
    int channel)
{
    Ray shadow_ray(start_point, normalize(end_point - start_point));

    float t_near, t_far;
    bool shade_vol = medium->intersect(shadow_ray, t_near, t_far);
    if (!shade_vol)
    {
        return 1;
    }

    float max_t = f_min(t_far, distance(start_point, end_point));

    int nsamp = 2;
    float count = 0; // the samples that reached the destination

    for (int n = 0; n < nsamp; n++)
    {
        float dist = t_near;

        for (;;)
        {
            dist += medium->sampleFreeDistance(rand01(), channel);
            if (dist >= max_t)
            {
                count += 1;
                break;
            }
            vec3 pos = shadow_ray.at(dist);

            if (medium->isNotNullCollision(pos, rand01()))
            {
                break;
            }
        }
    }
    return count / nsamp;
}

float volumetricPathTacing(
    const Ray& ray,
    int max_depth,
    int channel)
{
    if (0)
    {
        vec3 p;
        return light.Li(ray.o, ray.d, p)[channel];
    }

    Ray cr(ray);
    float radiance = (0.0f);
    float throughput = (1.0f);

    for (int depth = 0; depth < max_depth; depth++)
    {
        float t_near, t_far;
        bool shade_vol = medium->intersect(cr, t_near, t_far); // TODO : use scene to manage volume and light

        if (!shade_vol)
        {
            // using direct lighting otherwise
            // not attenuated since not in volume
            if (depth > 0 || P::show_bg)
            {
                vec3 lpos;
                radiance += (light.Li(cr.o, cr.d, lpos)[channel] * throughput);
            }
            break;
        }

        /// woodcock tracking / delta tracking
        vec3 pos = cr.at(t_near); // current position
        float dist = t_near;

        bool through = false;
        for (;;)
        {
            dist += medium->sampleFreeDistance(rand01(), channel);
            if (dist >= t_far)
            {
                through = true; // transmitted through the volume, probability is 1-exp(-optical_thickness)
                break;
            }
            pos = cr.at(dist);
            if (medium->isNotNullCollision(pos, rand01()))
            {
                break;
            }
        }

        // probability is exp(-optical_thickness)
        if (through)
        {
            if (depth > 0 || P::show_bg)
            {
                vec3 lpos;
                radiance += (light.Li(cr.o, cr.d, lpos)[channel] * throughput);
            }
            break;
        }

        // incoming radiance is scattered by sigma_s, and the line sampling pdf by delta tracking
        // is sigma_t * exp(-optical_thickness), medium attenuation is exp(-optical_thickness),
        // by cancelling only the albedo (= sigma_s / sigma_t) remains
        throughput *= medium->albedo(channel);

        Frame frame(cr.d);

        // MIS direct lighting
        if (1) // sample light
        {
            vec3 Li;
            float pdfW;
            vec3 lpos, ldir;
            light.sample(pos, pdfW, lpos, ldir, Li);

            float pdfW_phase = phase->evaluate(frame, ldir);
            float mis_weight = 1.0f;// misWeightBalanceHeuristic(pdfW, pdfW_phase);

            float transmittance = stochasticTransmittance(pos, lpos, channel);
            float attenuated_radiance = Li[channel] * transmittance;
            radiance += (attenuated_radiance * throughput) *
                (phase->evaluate(frame, ldir) / pdfW) * mis_weight;
        }
        if (0) // sample phase function
        {
            vec3 lpos;
            vec3 ldir = phase->sample(frame, rand01(), rand01());
            vec3 Li = light.Li(pos, ldir, lpos);
            float pdfW = phase->evaluate(frame, ldir);

            vec3 temp1, temp2;
            float pdfW_light = light.pdf(pos, ldir, temp1, temp2);
            float mis_weight = misWeightBalanceHeuristic(pdfW, pdfW_light);

            float transmittance = stochasticTransmittance(pos, lpos, channel);
            float attenuated_radiance = Li[channel] * transmittance;
            radiance += (attenuated_radiance * throughput) *
                mis_weight; // phase function and pdfW cancelled
        }

        // scattered direction
        vec3 new_dir = phase->sample(frame, rand01(), rand01());
        cr = Ray(pos, new_dir);
    }

    if (radiance != radiance)
    {
        radiance = 0.0f;
    }

    return radiance;
}

void render()
{
    FrameBuffer framebuffer(P::width, P::height);

    for (int k = 0; k < P::spp; k++)
    {
        for (int j = 0; j < P::height; j++)
        {
            printf("\rfinished %.2f%% - %d / %d     ", 100.0f * float(k + 1) / float(P::spp), j, P::height);

#pragma omp parallel for
            for (int i = 0; i < P::width; i++)
            {
                framebuffer.buffer()[i + j * P::width] += vec3(
                    volumetricPathTacing(cam->pixelRay(i, j), P::trace_depth, 0),
                    volumetricPathTacing(cam->pixelRay(i, j), P::trace_depth, 1),
                    volumetricPathTacing(cam->pixelRay(i, j), P::trace_depth, 2)
                    );
            }
        }

        FrameBuffer copy = framebuffer;
        copy.scaleBrightness(1.0f / (k + 1));
        copy.dump();
    }

    framebuffer.scaleBrightness(1.0f / P::spp);
    framebuffer.dump();

    printf("\rrendering finished\n");
}

int main()
{
    P::param();
    cam.reset(new Camera(P::camera_origin, P::camera_lookat, P::fov, P::width, P::height));

    printf("DEBUG:: begin program\n");

    vol.reset(new Texture3D(128, 1, vec3(-0.5)));
    vol->initChecker();
//     vol->initConstant(1.0f);
//     vol.load_binary("smoke.vol");

    medium.reset(new HeterogeneousMedium(vol, P::sigma_s, P::sigma_a));
    phase.reset(new HGPhaseFunction(P::HG_mean_cosine));

    printf("init complete\n");
    printf("rendering\n");
    render();

    return 0;
}