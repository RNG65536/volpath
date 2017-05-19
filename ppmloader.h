#ifndef ppmLoader_h__
#define ppmLoader_h__

vec3 envmap_sample_dir(vec3 r, const float *d_pixels,
    int width, int height)
{
    float xi = ((r.x >= 0 ? atanf(r.z / r.x) : atanf(r.z / r.x) + M_PI) + M_PI_2) / (2 * M_PI), yi = acosf(r.y) / M_PI;

    int x_ = int(xi*(width - 1) + .5);
    int y_ = int(yi*(height - 1) + .5);
    if (x_<0 || x_>width - 1 || y_<0 || y_>height - 1)
        return vec3(0, 0, 0);
    int offset = (x_ + y_*width) * 3;
    return vec3(d_pixels[offset], d_pixels[offset + 1], d_pixels[offset + 2]);
}


struct dir_light{
    vec3 direction;
    vec3 radiance;
    float pdf;
};

typedef enum{ CLAMP, CYCLIC, NONE } kind_;
class ppmLoader
{
private:
    std::vector<float> th_h_rawPixels;

public:
    char *filename;
    int width, height;
    float *p_d_rawPixels;
    float *p_h_rawPixels;

    float *pu_host, *pu_device;
    float *Pu_host, *Pu_device;
    float *pv_host, *pv_device;
    float *Pv_host, *Pv_device;

    ppmLoader() :p_d_rawPixels(NULL), p_h_rawPixels(NULL){}
    ~ppmLoader(){
//         if (p_d_rawPixels){
//             std::cout << "Freeing envmap" << std::endl;
//             memset((void*)p_d_rawPixels, 0, th_h_rawPixels.size()*sizeof(float));
//             free(p_d_rawPixels);
//         }

        delete[] pu_host;
        delete[] Pu_host;
        delete[] pv_host;
        delete[] Pv_host;
        free(pu_device);
        free(Pu_device);
        free(pv_device);
        free(Pv_device);
    }
    int openImageFile_hdr(const char *filename, float scale = 1.0f, float rotate_degree = 0.0f, int count = 10000);
    vec3 getRGB(float x, float y, kind_ k = NONE);
    vec3 getRGBdevice(float x, float y, kind_ k = NONE);
    void toDevice();
    void init_pdf();

    void presample_envmap(int count);
    std::vector<dir_light> h_light;
    std::vector<dir_light> d_light;
};

int ppmLoader::openImageFile_hdr(
    const char *filename,
    float scale,
    float rotate_degree,
    int count)
{
    width = 600, height = 300;
    std::ifstream ifs(filename, std::ios::in | std::ios::binary);
    if (ifs.good()) {
        ifs.read(reinterpret_cast<char*>(&width), sizeof(int));
        ifs.read(reinterpret_cast<char*>(&height), sizeof(int));
    }
    else {
        fprintf(stderr, "Can't open file %s\n", filename);
        ifs.close();
    }

    th_h_rawPixels.resize(width * height * 3);
    p_h_rawPixels = rawPtr(th_h_rawPixels);
    if (!p_h_rawPixels) return 1;

    if (ifs.good()) {
        ifs.read(reinterpret_cast<char*>(p_h_rawPixels), width*height * 3 * sizeof(float));
    }
    else {
        float *p = p_h_rawPixels;
        for (unsigned j = 0; j < height; ++j) {
            for (unsigned i = 0; i < width; ++i) {
                float hl = expf(-100 * (sq((i - width*0.0) / width) + sq((j - height*0.5) / height))) * 30 + 0.3;
                p[0] = hl*.98;
                p[1] = hl*.95;
                p[2] = hl*.92;
                p += 3;
            }
        }
    }

    std::vector<float> temp = th_h_rawPixels;
    int rshift(rotate_degree*width / 360.0f);
    printf("shift by %d pixels\n", rshift);
    for (int j = 0; j < height; j++){
        for (int i = 0; i < width; i++){
            int offset = i + j*width;
            int offset2 = (i - rshift + width) % width + j*width;
            th_h_rawPixels[offset * 3] = temp[offset2 * 3];
            th_h_rawPixels[offset * 3 + 1] = temp[offset2 * 3 + 1];
            th_h_rawPixels[offset * 3 + 2] = temp[offset2 * 3 + 2];
        }
    }

    ifs.close();

    for (int n = 0; n < th_h_rawPixels.size(); n++){
        th_h_rawPixels[n] *= scale;
    }

    init_pdf();
    presample_envmap(count);

    toDevice();
    return 0;
}

void ppmLoader::toDevice()
{
    size_t total = th_h_rawPixels.size()*sizeof(float);
    p_d_rawPixels = (float*)malloc(total);
    memcpy((void*)p_d_rawPixels, (void*)p_h_rawPixels, total);

    pu_device = (float*)malloc((height)*sizeof(float));
    Pu_device = (float*)malloc((height + 1)*sizeof(float));
    pv_device = (float*)malloc((height*width)*sizeof(float));
    Pv_device = (float*)malloc((height*(width + 1))*sizeof(float));
    memcpy((void*)pu_device, (void*)pu_host, (height)*sizeof(float));
    memcpy((void*)Pu_device, (void*)Pu_host, (height + 1)*sizeof(float));
    memcpy((void*)pv_device, (void*)pv_host, (height*width)*sizeof(float));
    memcpy((void*)Pv_device, (void*)Pv_host, (height*(width + 1))*sizeof(float));
}

template <typename T>
static inline T clamp(T x, T a, T b)
{
    return x<a ? a : x>b ? b : x;
}

vec3 ppmLoader::getRGB(float x, float y, kind_ k)
{
    int x_ = int(x*(width - 1) + .5);
    int y_ = int(y*(height - 1) + .5);
    int offset;

    switch (k){
    case CLAMP:
        x_ = clamp(x_, 0, width - 1);
        y_ = clamp(y_, 0, height - 1);
        offset = (x_ + y_*width) * 3;
        break;
    case CYCLIC:
        //         x_ %= width;
        //         y_ %= height;
        x_ = (x_ + width * 1000) % width;
        y_ = (y_ + height * 1000) % height;
        offset = (x_ + y_*width) * 3;
        break;
    case NONE:
    default:
        if (x_<0 || x_>width - 1 || y_<0 || y_>height - 1)
            return vec3(0, 0, 0);
        offset = (x_ + y_*width) * 3;
        break;
    }
    return vec3(p_h_rawPixels[offset],
        p_h_rawPixels[offset + 1], p_h_rawPixels[offset + 2]);
}

/************************************************************************/
/* Inversion method for importance sampling                             */
/************************************************************************/
static void precompute1D(float *f, float *pf, float *Pf, int nf){
    float I = 0.0f;
    for (int n = 0; n < nf; n++){
        I += f[n];
    }
    for (int n = 0; n < nf; n++){
        pf[n] = f[n] / I;
    }
    Pf[0] = 0.0f;
    for (int n = 1; n < nf; n++){
        Pf[n] = Pf[n - 1] + pf[n - 1];
    }
    Pf[nf] = 1.0f;
}

static void precompute2D(
    float *f, float *pu, float *Pu, float *pv, float *Pv,
    int nu, int nv)
{
    std::vector<float> rowsum(nu);
    for (int u = 0; u < nu; u++){
        //         precompute1D(f+u*nv, pv+u*nv, Pv+u*nv, nv);
        precompute1D(f + u*nv, pv + u*nv, Pv + u*(nv + 1), nv);//bug fixed
        float temp = 0.0f;
        for (int n = 0; n < nv; n++){
            temp += (f + u*nv)[n];
        }
        rowsum[u] = temp;
    }
    precompute1D(rawPtr(rowsum), pu, Pu, nu);
}

static void sample1D(float *pf, float *Pf, float unif, int nf,
    float &x, float &p)
{
    int i = 0;
    for (; i < nf; i++){
        if (Pf[i] <= unif && unif < Pf[i + 1] || i == nf - 1)//avoid overflow by additional check
            break;
    }
    float t = Pf[i + 1] > Pf[i] ?
        (Pf[i + 1] - unif) / (Pf[i + 1] - Pf[i])
        : 0;

    x = (1 - t) * i + t * (i + 1);
    p = pf[i];
}

static void sample2D(float *pu, float *Pu, float *pv, float *Pv,
    float unif1, float unif2, int nu, int nv,
    float &u, float &v, float &pdf)
{
    float pdfu, pdfv;
    sample1D(pu, Pu, unif1, nu, u, pdfu);
    sample1D(pv + int(u)*nv, Pv + int(u)*(nv + 1), unif2, nv, v, pdfv);
    pdf = pdfu * pdfv;
}

void ppmLoader::init_pdf(){
    std::vector<float> luminance(width*height);
    float luminance_sum = 0;
    for (int n = 0; n < width*height; n++){
        int offset = n * 3;
        luminance[n] =
            p_h_rawPixels[offset] * 0.2126
            + p_h_rawPixels[offset + 1] * 0.7152
            + p_h_rawPixels[offset + 2] * 0.0722;
        luminance_sum += luminance[n];
    }

    pu_host = new float[height];
    Pu_host = new float[height + 1];
    pv_host = new float[height*width];
    Pv_host = new float[height*(width + 1)];

    precompute2D(rawPtr(luminance), pu_host, Pu_host,
        pv_host, Pv_host, height, width);

    //test
    float z0 = 0;
    int test_count = 10000;
    float integral = 0;
    for (int n = 0; n < test_count; n++){
        float test_u, test_v, test_pdf;
        sample2D(pu_host, Pu_host, pv_host, Pv_host, rand01(), rand01(), height, width,
            test_u, test_v, test_pdf);
        z0 += test_pdf;
        integral += (luminance[int(test_u)*width + int(test_v)] / luminance_sum) / test_pdf;

    }
    integral /= test_count;//integrand is luminance[i]/luminance_sum, should integrate to 1
    printf("z0=%f\n", z0);
    printf("z1=%f\n", float(test_count) / (width*height));
    printf("int=%f\n", integral);
}

void ppmLoader::presample_envmap(int count)
{
    int N = int(sqrtf(count));
    int n_light = N*N;
    std::vector<float> e1, e2;
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            e1.push_back((i + rand01()) / float(N));
            e2.push_back((j + rand01()) / float(N));
        }
    }

    for (int n = 0; n < n_light; n++)
    {
        float test_u, test_v, test_pdf;
        sample2D(pu_host, Pu_host, pv_host, Pv_host, e1[n], e2[n], height, width,
            test_u, test_v, test_pdf);

        int pixel_idx = int(test_u)*width + int(test_v);
        vec3 sample_radiance(
            p_h_rawPixels[pixel_idx * 3],
            p_h_rawPixels[pixel_idx * 3 + 1],
            p_h_rawPixels[pixel_idx * 3 + 2]
            );

        float theta = (int(test_u) + 0.5) / float(height) * M_PI;
        float phi = (int(test_v) + 0.5) / float(width) * (2.0f * M_PI);//azimuthal
        vec3 sample_dir = vec3(sinf(theta)*sinf(phi), cosf(theta), sinf(theta)*-cosf(phi));

        float sample_pdf = test_pdf * ((float(width)*float(height)) / (2.0f*M_PI*M_PI*sinf(theta)));//warpping
        dir_light dl;
        dl.direction = sample_dir;
        dl.radiance = sample_radiance;
        dl.pdf = sample_pdf;

        h_light.push_back(dl);
    }

    d_light = h_light;

    printf(".#light=..%d...\n", d_light.size());
}

ppmLoader ppm;
#endif // ppmLoader_h__
