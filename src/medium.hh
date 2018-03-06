#ifndef medium_h__
#define medium_h__

class HGPhaseFunction
{
    float g;

    // perfect inversion, pdf matches evaluation exactly
    vec3 sample(float rnd0, float rnd1) const
    {
        float cos_theta;
        if (std::abs(g) > 1e-6f)
        {
            float s = 2.0f * rnd0 - 1.0f;
            float f = (1.0f - g * g) / (1.0f + g * s);
            cos_theta = (0.5f / g) * (1.0f + g * g - f * f);
            cos_theta = f_max(0.0f, f_min(1.0f, cos_theta));
        }
        else
        {
            cos_theta = 2.0f * rnd0 - 1.0f;
        }
        float sin_theta = std::sqrt(1.0f - cos_theta * cos_theta);
        float phi = 2.0f * M_PI * rnd1;
        vec3 ret = (vec3(std::cos(phi) * sin_theta, std::sin(phi) * sin_theta, cos_theta));
        return ret;
    }

    float evaluate(float cos_theta) const
    {
        return (1.0f - g * g) / (4.0f * M_PI * std::pow(1.0f + g * g - 2 * g * cos_theta, 1.5f));
    }

public:
    HGPhaseFunction(float g)
        : g(g)
    {

    }

    vec3 sample(const Frame& frame, float rnd0, float rnd1) const
    {
        vec3 s = sample(rnd0, rnd1);
        return frame.toWorld(s);
    }

    float evaluate(const Frame& frame, const vec3& dir) const
    {
        float cos_theta = dot(frame.normal(), dir);
        return evaluate(cos_theta);
    }
};

class HeterogeneousMedium
{
    std::unique_ptr<Texture3D> m_density;
    vec3 m_sigma_s, m_sigma_a, m_sigma_t;
    vec3 m_albedo;

    float m_max_density;
    float m_inv_max_density;
    vec3 m_max_sigma_t; // sigma is directly modulated by density
    vec3 m_inv_max_sigma_t; // sigma is directly modulated by density
    float m_max_sigma_t_max_compoennt; // usually as majorant
    float m_inv_max_sigma_t_max_component;

public:
    HeterogeneousMedium(std::unique_ptr<Texture3D>& vol, const vec3& sigma_s, const vec3& sigma_a)
    {
        m_density = std::move(vol);
        float max_density = m_density->maxValue();
        printf("maximum density is %f\n", max_density);
        m_max_density = max_density;
        m_inv_max_density = 1.0f / m_max_density;

        m_sigma_s = sigma_s;
        m_sigma_a = sigma_a;

        m_sigma_t = m_sigma_s + m_sigma_a;
        m_albedo = m_sigma_s / m_sigma_t;
        m_max_sigma_t = m_sigma_t * max_density;
        m_inv_max_sigma_t = vec3(1.0f) / m_max_sigma_t;

        m_max_sigma_t_max_compoennt =
            f_max(f_max(m_max_sigma_t.x, m_max_sigma_t.y), m_max_sigma_t.z);
        m_inv_max_sigma_t_max_component = 1.0f / m_max_sigma_t_max_compoennt;
    }

    vec3 sigmaT(const vec3& pos) const
    {
        return m_sigma_t * m_density->fetch(pos);
    }
    float sigmaT(const vec3& pos, int channel) const
    {
        return m_sigma_t[channel] * m_density->fetch(pos);
    }

    vec3 sigmaS(const vec3& pos) const
    {
        return m_sigma_s * m_density->fetch(pos);
    }
    float sigmaS(const vec3& pos, int channel) const
    {
        return m_sigma_s[channel] * m_density->fetch(pos);
    }

    vec3 sigmaA(const vec3& pos) const
    {
        return m_sigma_a * m_density->fetch(pos);
    }
    float sigmaA(const vec3& pos, int channel) const
    {
        return m_sigma_a[channel] * m_density->fetch(pos);
    }

    vec3 albedo() const
    {
        return m_albedo;
    }
    float albedo(int channel) const
    {
        return m_albedo[channel];
    }

    float density(const vec3& pos) const
    {
        return m_density->fetch(pos);
    }

    //////////////////////////////////////////////////////////////////////////

    float maxDensity() const
    {
        return m_max_density;
    }
    float invMaxDensity() const
    {
        return m_inv_max_density;
    }
    float maxSigmaT(int channel) const
    {
        return m_max_sigma_t[channel];
    }
    float invMaxSigmaT(int channel) const
    {
        return m_inv_max_sigma_t[channel];
    }
    float maxSigmaTMaxComponent() const
    {
        return m_max_sigma_t_max_compoennt;
    }
    float invMaxSigmaTMaxComponent() const
    {
        return m_inv_max_sigma_t_max_component;
    }

    float sampleFreePath(float rnd, float inv_sigma) const
    {
        return -std::log(rnd + FLT_MIN) * inv_sigma; // log(0) -> -INF
    }
    //////////////////////////////////////////////////////////////////////////
//     // continue only if null collision
//     bool isNotNullCollision(const vec3& pos, float rnd) const
//     {
//         // sigma_t / sigma_t_max
//         return rnd < m_density->fetch(pos) * m_inv_max_density;
//     }
//     float sampleFreeDistance(float rnd, int channel) const
//     {
// //         return -std::log(1.0f - rnd) * m_inv_sigma_t_max[channel];
//         return -std::log(rnd + FLT_MIN) * m_inv_max_sigma_t[channel]; // log(0) -> -INF
//     }

    bool intersect(const Ray& r, float& t_near, float& t_far) const
    {
        return intersect_vol(r, m_density->min, m_density->max, t_near, t_far);
    }
};

#endif // medium_h__
