// assuming the same sigma and albedo for all the channels
class VolpathCloud : public VolumetricPathTracer
{
    float deltaTrackingTransmittanceEstimatorChannel(
        const vec3& start_point,
        const vec3& end_point,
        int channel)
    {
        Ray shadow_ray(start_point, normalize(end_point - start_point));

        float t_near, t_far;
        bool shade_vol = medium->intersect(shadow_ray, t_near, t_far);
        if (!shade_vol)
        {
            return 1.0f;
        }

        float inv_sigma = 1.0f / medium->maxSigmaT(channel);
        float max_t = f_min(t_far, distance(start_point, end_point));

        float dist = t_near;
        for (;;)
        {
            dist += medium->sampleFreePath(rand01(), inv_sigma);
            if (dist >= max_t)
            {
                break;
            }
            vec3 pos = shadow_ray.at(dist);

            if (rand01() < medium->sigmaT(pos)[channel] * inv_sigma)
            {
                break;
            }
        }
        return float(dist >= max_t);
    }

public:
    vec3 trace(
        const Ray& ray,
        int max_depth,
        Stats& stats)
    {
        if (0)
        {
            vec3 p;
            return light.Li(ray.o, ray.d, p);
        }

        Ray cr(ray);
        vec3 radiance = vec3(0.0f);
        float throughput = 1.0f;

        int depth;
        for (depth = 0; depth < max_depth; depth++)
        {
//             if (1)
//             {
//                 const float RR_continue_prob = depth < 100 ? 1.0f : 0.5f;
//                 if (rand01() > RR_continue_prob)
//                 {
//                     break;
//                 }
//                 throughput /= RR_continue_prob;
//             }

            float t_near, t_far;
            bool shade_vol = medium->intersect(cr, t_near, t_far); // TODO : use scene to manage volume and light

            // this has nothing to do with the medium
            if (!shade_vol)
            {
                if (depth > 0 || P::show_bg)
                {
                    vec3 lpos;
                    radiance += light.Li(cr.o, cr.d, lpos) * throughput;
                }

                break;
            }

            /// woodcock tracking / delta tracking
            vec3 pos = cr.at(t_near); // current position
            float dist = t_near;
            int sampling_channel = 0;
            float inv_sigma = medium->invMaxSigmaT(sampling_channel);

            bool through = false;
            // delta tracking scattering event sampling
            for (;;)
            {
                dist += medium->sampleFreePath(rand01(), inv_sigma);
                pos = cr.at(dist);
                if (dist >= t_far)
                {
                    through = true; // transmitted through the volume, probability is 1-exp(-optical_thickness)
                    break;
                }
                if (rand01() < medium->sigmaT(pos)[sampling_channel] * inv_sigma)
                {
                    break;
                }
            }

            // probability is exp(-optical_thickness)
            if (through)
            {
                if (depth > 0 || P::show_bg) // if direct lighting does not consider the background, seems lowest variance
                {
                    vec3 lpos;
                    bool show_sun = (depth == 0);
                    radiance += light.Li(cr.o, cr.d, lpos, show_sun) * throughput;
                }

                break;
            }

            // incoming radiance is scattered by sigma_s, and the line sampling pdf by delta tracking
            // is sigma_t * exp(-optical_thickness), medium attenuation is exp(-optical_thickness),
            // by cancelling only the albedo (= sigma_s / sigma_t) remains
            throughput *= medium->albedo()[sampling_channel];

            Frame frame(cr.d);

            // direct lighting
            if (1) // sample light
            {
                vec3 Li;
                float pdfW;
                vec3 lpos, ldir;
                light.sample(pos, pdfW, lpos, ldir, Li);

                float transmittance = deltaTrackingTransmittanceEstimatorChannel(pos, lpos, 0);
                if (transmittance > 0)
                {
                    vec3 attenuated_radiance = Li * transmittance;
                    radiance += (attenuated_radiance * throughput) * (phase->evaluate(frame, ldir) / pdfW);
                }
            }

            // scattered direction
            vec3 new_dir = phase->sample(frame, rand01(), rand01());
            cr = Ray(pos, new_dir);
        }


        stats.min_depth = std::min(stats.min_depth, depth);
        stats.max_depth = std::max(stats.max_depth, depth);
        stats.depth_sum += depth;

        if (radiance != radiance)
        {
            radiance = vec3(0.0f);
        }

        return radiance;
    }
};
