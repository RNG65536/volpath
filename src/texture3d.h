#ifndef texture3d_h__
#define texture3d_h__

inline int index_xyz(int x, int y, int z, int N, int NN) { return x + y * N + z * NN; }
inline int index_xzy(int x, int y, int z, int N, int NN) { return x + z * N + y * NN; }
inline int index_yxz(int x, int y, int z, int N, int NN) { return y + x * N + z * NN; }
inline int index_yzx(int x, int y, int z, int N, int NN) { return y + z * N + x * NN; }
inline int index_zxy(int x, int y, int z, int N, int NN) { return z + x * N + y * NN; }
inline int index_zyx(int x, int y, int z, int N, int NN) { return z + y * N + x * NN; }
#define index_convention index_xyz

class Texture3D
{
    void scaleIntensity(float s)
    {
        for (auto& d : data)
        {
            d *= s;
        }
    }

public:
    int N, NN, total;
    float length_x;
    vec3 min, max;
    std::vector<float> data;

    Texture3D()
        : length_x(1), min(vec3(-0.5)), max(min + length_x)
    {
        resize(1, 1, 1);
    }
    Texture3D(int n, float length_x = 1, const vec3& base_pos = vec3(-0.5))
        : length_x(length_x), min(base_pos), max(min + length_x)
    {
        resize(n, n, n);
    }
    void resize(int x, int y, int z)
    {
        N = x;
        NN = N * N;
        total = N * NN;
        data.resize(total);
        data.assign(total, float(0.0));
    }
    void loadBinary(const char *fn)
    {
        struct VolumeHeader
        {
            uint8_t pad0[8];
            int mx, my, mz;
            int channel;
            uint8_t pad1[24];
        }; // totally 48 bytes
        VolumeHeader header;

        FILE *fp;
        fopen_s(&fp, fn, "rb");
        fread(&header, 1, sizeof(VolumeHeader), fp);
        resize(header.mx, header.my, header.mz);
        fread(&data[0], sizeof(float), total, fp);
        fclose(fp);
        printf("loaded %s <%d,%d,%d>\n", fn, header.mx, header.my, header.mz);
    }
    void initChecker()
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                for (int k = 0; k < N; k++)
                {
                    int xi = (int)ceil(5.0*(i + 1) / N);
                    int yi = (int)ceil(5.0*(j + 1) / N);
                    int zi = (int)ceil(5.0*(k + 1) / N);
                    data[index(i, j, k)] = float((xi + yi + zi) & 0x01);
                }
            }
        }
    }
    void initCheckerSmooth()
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                for (int k = 0; k < N; k++)
                {
                    float xi = std::sin(1.0 * (i + 1) / N * M_PI);
                    float yi = std::sin(1.0 * (j + 1) / N * M_PI);
                    float zi = std::sin(1.0 * (k + 1) / N * M_PI);
                    data[index(i, j, k)] = f_max(0.0f, xi * yi * zi);
                }
            }
        }
    }
    void initConstant(float f)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                for (int k = 0; k < N; k++)
                {
                    data[index(i, j, k)] = f;
                }
            }
        }
    }

    float& operator[] (int n) { return data[n]; }
    const float& operator[] (int n) const { return data[n]; }
    float& operator() (int i, int j, int k) { return data[index(i, j, k)]; }
    const float& operator() (int i, int j, int k) const { return data[index(i, j, k)]; }

    bool isOutside(const vec3& pos) const
    {
        return pos.x < min.x || pos.y < min.y || pos.z < min.z
            || pos.x > max.x || pos.y > max.y || pos.z > max.z;
    }

    vec3 toWorld(int i, int j, int k)
    {
        return vec3((i + 0.5f) / (float)N,
                    (j + 0.5f) / (float)N,
                    (k + 0.5f) / (float)N) * length_x + min;
    }

    vec3 toLocal(const vec3& pos) const
    {
        return (pos - min) / length_x;
    }

    inline int index(int x, int y, int z) const
    {
        return index_convention(x, y, z, N, NN);
    }

    void binarize(float ratio)
    {
        float _min = *std::min_element(data.begin(), data.end());
        float _max = *std::max_element(data.begin(), data.end());
        const float threshold = ratio * (_max - _min) + _min;
        std::transform(data.begin(), data.end(), data.begin(), [&threshold](float x)->float
        {
            return x < threshold ? 0.0f : 1.0f;
        });
    }

    void unitize()
    {
        auto result = std::minmax_element(data.begin(), data.end());
        float _min = *result.first;
        float _max = *result.second;
        printf("before is <%f, %f>\n", _min, _max);
        float inv = 1.0f / (_max - _min);
        std::transform(data.begin(), data.end(), data.begin(), [&inv, &_min](float x)->float
        {
            return  (x - _min) * inv;
        });
        _min = *std::min_element(data.begin(), data.end());
        _max = *std::max_element(data.begin(), data.end());
        printf("after is <%f, %f>\n", _min, _max);
    }

    float maxValue()
    {
        return *std::max_element(data.begin(), data.end());
    }

    float fetch(const vec3& pos) const
    {
        return 1.0f;
        if (isOutside(pos))
        {
            return 0;
        }

        vec3 p = toLocal(pos);
        float x_ = p.x;
        float y_ = p.y;
        float z_ = p.z;

        //trilinear
        x_ = x_ * N - 0.5f;
        y_ = y_ * N - 0.5f;
        z_ = z_ * N - 0.5f;

        int i_0 = (int)x_;
        int j_0 = (int)y_;
        int k_0 = (int)z_;
        int i_1 = i_0 + 1;
        int j_1 = j_0 + 1;
        int k_1 = k_0 + 1;

        //fixed
        i_0 = clampi(i_0, 0, N - 1);
        j_0 = clampi(j_0, 0, N - 1);
        k_0 = clampi(k_0, 0, N - 1);
        i_1 = clampi(i_1, 0, N - 1);
        j_1 = clampi(j_1, 0, N - 1);
        k_1 = clampi(k_1, 0, N - 1);

        float u1 = x_ - i_0;
        float v1 = y_ - j_0;
        float w1 = z_ - k_0;
        float u0 = 1 - u1;
        float v0 = 1 - v1;
        float w0 = 1 - w1;

        return u0 * (v0 * (w0 * (data[index(i_0, j_0, k_0)])
            +              w1 * (data[index(i_0, j_0, k_1)]))
            +        v1 * (w0 * (data[index(i_0, j_1, k_0)])
            +              w1 * (data[index(i_0, j_1, k_1)])))
            +  u1 * (v0 * (w0 * (data[index(i_1, j_0, k_0)])
            +              w1 * (data[index(i_1, j_0, k_1)]))
            +        v1 * (w0 * (data[index(i_1, j_1, k_0)])
            +              w1 * (data[index(i_1, j_1, k_1)])))
            ;
    }
};

#endif // texture3d_h__
