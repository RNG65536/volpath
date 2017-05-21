class Volpath : public VolumetricPathTracer
{
public:
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

            if (rand01() <
                ( medium->sigmaT(pos)[channel] / medium->maxSigmaT(channel) )
                )
            {
                break;
            }
        }
        return float(dist >= max_t);
    }

    float deltaTrackingEstimatorCompensationChannel(
        const vec3& start_point,
        const vec3& end_point,
        int sampling_channel,
        int channel)
    {
        Ray shadow_ray(start_point, normalize(end_point - start_point));

        float t_near, t_far;
        bool shade_vol = medium->intersect(shadow_ray, t_near, t_far);
        if (!shade_vol)
        {
            return 1.0f;
        }

        float inv_sigma = 1.0f / medium->maxSigmaT(sampling_channel);
        float max_t = f_min(t_far, distance(start_point, end_point));
        float Tr = 1.0f;

        float dist = t_near;
        for (;;)
        {
            dist += medium->sampleFreePath(rand01(), inv_sigma);
            if (dist >= max_t)
            {
                break;
            }
            vec3 pos = shadow_ray.at(dist);

            Tr *= 1.0f - (medium->sigmaT(pos)[channel] - medium->sigmaT(pos)[sampling_channel]) / (medium->maxSigmaT(channel));
//             if (rand01() < (medium->sigmaT(pos)[channel] - medium->sigmaT(pos)[sampling_channel]) / (medium->maxSigmaT(channel)))
//             {
//                 break;
//             }
        }
//         return float(dist >= max_t);
        return Tr;
    }

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
        vec3 throughput = vec3(1.0f);

        int depth;
        for (depth = 0; depth < max_depth; depth++)
        {
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
            int sampling_channel = std::min(int(rand01() * 3.0f), 2);
            float sigma = medium->maxSigmaT(sampling_channel);
            float inv_sigma = 1.0 / sigma;
            float sample_dist;

            // sampling exponentially
            bool through = false;
            for (;;)
            {
                sample_dist = medium->sampleFreePath(rand01(), inv_sigma);
                dist += sample_dist;
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

            vec3 Tr;
            for (int channel = 0; channel < 3; channel++)
            {
                Tr[channel] = deltaTrackingEstimatorCompensationChannel(cr.at(t_near), pos, sampling_channel, channel);
            }
            throughput *= Tr * (medium->sigmaT(pos) / medium->sigmaT(pos)[sampling_channel]);


            // incoming radiance is scattered by sigma_s, and the line sampling pdf by delta tracking
            // is sigma_t * exp(-optical_thickness), medium attenuation is exp(-optical_thickness),
            // by cancelling only the albedo (= sigma_s / sigma_t) remains
            throughput *= medium->albedo();

            Frame frame(cr.d);

            // MIS direct lighting
            if (1) // sample light
            {
                vec3 Li;
                float pdfW;
                vec3 lpos, ldir;
                light.sample(pos, pdfW, lpos, ldir, Li);
                vec3 attenuated_radiance;

                for (int channel = 0; channel < 3; channel++)
                {
                    float transmittance = deltaTrackingEstimatorChannel(pos, lpos, channel);
                    attenuated_radiance[channel] = Li[channel] * transmittance;
                }
                radiance += (attenuated_radiance * throughput)
                    * (phase->evaluate(frame, ldir) / pdfW);
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
