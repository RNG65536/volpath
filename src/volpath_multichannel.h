class VolpathMultichannel : public VolumetricPathTracer
{
public:
    // stochastic transmittance
    float deltaTrackingEstimatorChannel(
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
        float inv_sigma = medium->invMaxSigmaT(channel);

        int nsamp = 2; // limited improvement
        float count = 0; // the samples that reached the destination

        for (int n = 0; n < nsamp; n++)
        {
            float dist = t_near;

            for (;;)
            {
                dist += medium->sampleFreePath(rand01(), inv_sigma);
                if (dist >= max_t)
                {
                    count += 1;
                    break;
                }
                vec3 pos = shadow_ray.at(dist);

                if (rand01() < medium->density(pos) * medium->invMaxDensity())
                {
                    break;
                }
            }
        }
        return count / nsamp;
    }

    float traceChannel(
        const Ray& ray,
        int max_depth,
        int channel,
        Stats& stats)
    {
        if (0)
        {
            vec3 p;
            return light.Li(ray.o, ray.d, p)[channel];
        }

        Ray cr(ray);
        float radiance = (0.0f);
        float throughput = (1.0f);

        int depth;
        for (depth = 0; depth < max_depth; depth++)
        {
            //         if (0)
            //         {
            //             const float RR_continue_prob = depth < 100 ? 1.0f : medium->albedo(channel);
            //             if (rand01() > RR_continue_prob)
            //             {
            //                 break;
            //             }
            //             throughput /= RR_continue_prob;
            //         }

            float t_near, t_far;
            bool shade_vol = medium->intersect(cr, t_near, t_far); // TODO : use scene to manage volume and light

            // this has nothing to do with the medium
            if (!shade_vol)
            {
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
            float inv_sigma = medium->invMaxSigmaT(channel);

            bool through = false;
            for (;;)
            {
                dist += medium->sampleFreePath(rand01(), inv_sigma);
                if (dist >= t_far)
                {
                    through = true; // transmitted through the volume, probability is 1-exp(-optical_thickness)
                    break;
                }
                pos = cr.at(dist);
                if (rand01() < medium->density(pos) * medium->invMaxDensity())
                {
                    break;
                }
            }

#define PASSIVE_BACKGROUND_SAMPLING 1

            // probability is exp(-optical_thickness)
            if (through)
            {
#if PASSIVE_BACKGROUND_SAMPLING
                if (depth > 0 || P::show_bg) // if direct lighting does not consider the background, seems lowest variance
#else
                if (depth == 0) // if using direct lighting from the background
#endif
                {
                    vec3 lpos;
                    bool show_sun = (depth == 0);
                    radiance += light.Li(cr.o, cr.d, lpos, show_sun)[channel] * throughput;
                }

                break;
            }

            // incoming radiance is scattered by sigma_s, and the line sampling pdf by delta tracking
            // is sigma_t * exp(-optical_thickness), medium attenuation is exp(-optical_thickness),
            // by cancelling only the albedo (= sigma_s / sigma_t) remains
            throughput *= medium->albedo(channel);

            Frame frame(cr.d);

            // MIS direct lighting
            float rad[2];
            float weight[2];
            if (1) // sample light
            {
                vec3 Li;
                float pdfW;
                vec3 lpos, ldir;
                light.sample(pos, pdfW, lpos, ldir, Li);

                float pdfW_phase = phase->evaluate(frame, ldir);
                float mis_weight = misWeightBalanceHeuristic(pdfW, pdfW_phase);

                float transmittance = deltaTrackingEstimatorChannel(pos, lpos, channel);
                float attenuated_radiance = Li[channel] * transmittance;
                rad[0] = (attenuated_radiance * throughput) * (phase->evaluate(frame, ldir) / pdfW);
                weight[0] = mis_weight;
            }
            if (0 &&
                !PASSIVE_BACKGROUND_SAMPLING) // sample phase function (allows background, but variance is higher than pure passive sampling)
            {
                vec3 lpos;
                vec3 ldir = phase->sample(frame, rand01(), rand01());
                vec3 Li = light.Li(pos, ldir, lpos);
                float pdfW = phase->evaluate(frame, ldir);

                vec3 temp1, temp2;
                float pdfW_light = light.pdf(pos, ldir, temp1, temp2);
                float mis_weight = misWeightBalanceHeuristic(pdfW, pdfW_light);

                float transmittance = deltaTrackingEstimatorChannel(pos, lpos, channel);
                float attenuated_radiance = Li[channel] * transmittance;
                rad[1] = (attenuated_radiance * throughput); // phase function and pdfW got cancelled
                weight[1] = mis_weight;
            }
            radiance += rad[0];// *weight[0];
            //         radiance += rad[1] * weight[1];

            // scattered direction
            vec3 new_dir = phase->sample(frame, rand01(), rand01());
            cr = Ray(pos, new_dir);
        }


        stats.min_depth = std::min(stats.min_depth, depth);
        stats.max_depth = std::max(stats.max_depth, depth);
        stats.depth_sum += depth;

        if (radiance != radiance)
        {
            radiance = 0.0f;
        }

        return radiance;
    }

    vec3 trace(const Ray& ray, int max_depth, Stats& stats)
    {
        float rad[3];
        rad[0] = traceChannel(ray, max_depth, 0, stats);
        rad[1] = traceChannel(ray, max_depth, 1, stats);
        rad[2] = traceChannel(ray, max_depth, 2, stats);
        return vec3(rad[0], rad[1], rad[2]);
    }
};
