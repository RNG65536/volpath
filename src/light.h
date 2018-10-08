#ifndef light_h__
#define light_h__

// TODO : fix light_pos :
// l = d * cos_theta_sample - sqrt(sq(r) - sq(d * sin_theta_sample))
// light_pos = position + light_dir * l

class Light
{
    const float distance = 1000.0f;
    const float brightness = 0.02f;

    // sphere light, mimics sun
    vec3 center = vec3(2, 2, 1);
//     vec3 center = vec3(-2, -2, -4);
    float radius = distance * 0.02f;
    vec3 emission = vec3(1.0f, 0.9f, 0.75f) * (0.8f * (distance * distance) / (radius * radius)); // sun radiance
//     vec3 emission = vec3(0.0f);
    Skydome sky;

public:
    Light() : sky(
        4.0f, // 3.0f, // 5.0f, // turbidity
        vec3(0.157f), // vec3(0.2f), // ground albedo
        M_PI * 0.1f, // sun elevation
        M_PI * (-0.2f + 0.5f), // sun azimuth
        0.02f) // brightness
    {
        center = sky.sunDirection() * distance;
    }
    void sample(vec3 position, float& pdfW, vec3& light_pos, vec3& light_dir, vec3& Li) const
    {
        Frame frame(center - position);

        float d2 = dot(position - center, position - center);
        float r2 = radius * radius;
        float l2 = d2 - r2;
        float cos_a_max = std::sqrt(l2 / d2);

        float rnd0 = rand01(), rnd1 = rand01();
        float cos_a = 1 - rnd0 + rnd0 * cos_a_max; // uniformly sample the cosine
        float sin_a = std::sqrt(1 - cos_a * cos_a);
        float phi = 2 * M_PI * rnd1;
        light_dir = frame.toWorld(vec3(std::cos(phi) * sin_a, std::sin(phi) * sin_a, cos_a));
        light_dir = normalize(light_dir);

        pdfW = 1.0f / (2.0f * M_PI * (1.0f - cos_a_max));
        Li = vec3(emission);
        float d = std::sqrt(d2);
        float dist_to_hit = d * cos_a - std::sqrt(r2 - sq(d * sin_a));
        light_pos = position + light_dir * dist_to_hit;
    }
    vec3 Li(const vec3& position, const vec3& direction, vec3& light_pos, bool include_sun = true) const
    {
        light_pos = position + direction * M_INFINITY;
        float d2 = dot(position - center, position - center);
        float l2 = d2 - radius * radius;
        if (l2 <= M_EPS) // inside sphere
        {
            return vec3(0.0f);
        }
        float cos_a_max = std::sqrt(l2 / d2);
        float cos_a = dot(direction, normalize(center - position));
        if (cos_a < cos_a_max)
        {
//             return vec3(0.0f);
            return sky.radiance(direction);
        }
        light_pos = position + direction * std::sqrt(l2);
        if (!include_sun) // because Li is used for both passive and aggressive lighting
        {
            return vec3(0.0f);
        }
        return vec3(emission); // todo : add intersection
    }
    float pdf(vec3& position, const vec3& direction, vec3& light_pos, vec3& Li) const
    {
        Li = vec3(0.0f);
        float d2 = dot(position - center, position - center);
        float l2 = d2 - radius * radius;
        if (l2 < 0) // inside sphere
        {
            return 0;
        }
        float cos_a_max = std::sqrt(l2 / d2);
        float cos_a = dot(direction, normalize(center - position));
        if (cos_a < cos_a_max)
        {
            return 0;
        }
        light_pos = position + direction * std::sqrt(l2);
        Li = vec3(emission);
        return 1.0f / (2.0f * M_PI * (1.0f - cos_a_max));
    }
};


#endif // light_h__
