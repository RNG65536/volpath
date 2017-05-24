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

#include "timer.h"
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
    int max_depth;

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

    width = 200;
    height = 200;
    spp = 1000;
    max_depth = 20000; // 1000; // 200; // 

    // sigma_t_prime = sigma_a + sigma_s * (1 - g)
    // for scattering dominated media, should scale with 1 / (1 - g) to approximate appearance
    density_scale = 500.0f; // 100.0f; //  40.0f;// 20.0f; // 5.0f; //  
    HG_mean_cosine = 0.877f;// 0.0f; //  0.7f;// 0.99f;// 

//     sigma_s = vec3(0.70f, 1.22f, 1.90f) * density_scale;
//     sigma_a = vec3(0.0014f, 0.0025f, 0.0142f) * density_scale;
//     sigma_a = vec3(0.0001f, 0.0001f, 0.0001f) * density_scale;

    sigma_s = vec3(1.0f, 1.0f, 1.0f) * density_scale;
    sigma_a = vec3(0.0f, 0.0f, 0.0f) * density_scale;

//     sigma_a = vec3(0.2f, 0.2f, 0.2f) * density_scale;
//     sigma_s = vec3(0.8f, 0.8f, 0.8f) * density_scale;

    // for compensation comparison
//     sigma_a = vec3(0.1f, 0.3f, 0.5f) * density_scale;
//     sigma_s = vec3(0.9f, 0.8f, 0.8f) * density_scale;
}

Light light;
std::unique_ptr<Camera> cam;
std::unique_ptr<HGPhaseFunction> phase;
std::unique_ptr<HeterogeneousMedium> medium;
std::unique_ptr<Texture3D> vol;

class VolumetricPathTracer
{
protected:
    class Stats
    {
    public:
        int depth_sum;
        int min_depth;
        int max_depth;

        void clear()
        {
            depth_sum = 0;
            min_depth = INT_MAX;
            max_depth = -INT_MAX;
        }

        void print() const
        {
            float mean_depth = (float)depth_sum / (float(P::width) * float(P::height) * 3.0f);
            cout << " depth mean, min, max : " << mean_depth << ", " << min_depth << ", " << max_depth << endl;
        }
    };

    Stats m_stats;

public:
    virtual vec3 trace(const Ray& ray, int max_depth, Stats& stats) = 0;

    void render()
    {
        FrameBuffer framebuffer(P::width, P::height);
        Timer timer;

        for (int k = 0; k < P::spp; k++)
        {
            m_stats.clear();

            for (int j = 0; j < P::height; j++)
            {
                printf("\rfinished %.2f%% - %d / %d     ", 100.0f * float(k + 1) / float(P::spp), j, P::height);

// #pragma omp parallel for
                for (int i = 0; i < P::width; i++)
                {
                    vec3 rad = trace(cam->pixelRay(i, j), P::max_depth, m_stats);
                    framebuffer.buffer()[i + j * P::width] += rad;
                }
            }

            m_stats.print();
            cout << " elapsed " << timer.elapsed() << " s, ";

            FrameBuffer copy = framebuffer;
            copy.scaleBrightness(1.0f / (k + 1));
            copy.dump();
        }

        framebuffer.scaleBrightness(1.0f / P::spp);
        framebuffer.dump();

        printf("\rrendering finished\n");
    }
};

#include "volpath_multichannel.h"
#include "volpath.h"
#include "volpath_cloud.h"

int main()
{
//     VolpathMultichannel volpath;
//     Volpath volpath;
    VolpathCloud volpath;

    P::param();
    cam.reset(new Camera(P::camera_origin, P::camera_lookat, P::fov, P::width, P::height));

    printf("DEBUG:: begin program\n");

    vol.reset(new Texture3D(128, 1, vec3(-0.5)));
    vol->initChecker();
//     vol->initCheckerSmooth();
//     vol->initConstant(1.0f);
//     vol.load_binary("smoke.vol");

    medium.reset(new HeterogeneousMedium(vol, P::sigma_s, P::sigma_a));
    phase.reset(new HGPhaseFunction(P::HG_mean_cosine));

    printf("init complete\n");
    printf("rendering\n");
    volpath.render();

    return 0;
}