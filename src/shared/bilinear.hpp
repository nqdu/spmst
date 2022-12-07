void bilinear(const float* restrict x, const float* restrict y,int nx,int ny,
                float x0,float y0,int &ix,int &iy,float* restrict coef);

float interp2d(const float* restrict x, const float* restrict y,
        const float* restrict z,int nx,int ny,float x0,float y0);