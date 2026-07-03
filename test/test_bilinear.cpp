#include "shared/bilinear.hpp"
#include <gtest/gtest.h>

TEST(Bilinear, ExactAtGridNode)
{
    float x[3] = {0.0, 1.0, 2.0};
    float y[3] = {0.0, 1.0, 2.0};
    float z[9] = {
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
        7.0, 8.0, 9.0,
    };

    // interpolation at a grid node should reproduce the exact value
    EXPECT_NEAR(interp2d(x,y,z,3,3,1.0,1.0), 5.0, 1e-5);
    EXPECT_NEAR(interp2d(x,y,z,3,3,0.0,0.0), 1.0, 1e-5);
    EXPECT_NEAR(interp2d(x,y,z,3,3,2.0,2.0), 9.0, 1e-5);
}

TEST(Bilinear, MidpointAverage)
{
    float x[2] = {0.0, 1.0};
    float y[2] = {0.0, 1.0};
    float z[4] = {
        0.0, 10.0,
        20.0, 30.0,
    };

    // midpoint of a bilinear patch is the average of the 4 corners
    EXPECT_NEAR(interp2d(x,y,z,2,2,0.5,0.5), 15.0, 1e-4);
}

TEST(Bilinear, ClampsOutsideDomain)
{
    float x[3] = {0.0, 1.0, 2.0};
    float y[3] = {0.0, 1.0, 2.0};

    int ix,iy;
    float coef[4];

    // requesting a point left/below the domain should clamp to the first cell
    bilinear(x,y,3,3,-5.0,-5.0,ix,iy,coef);
    EXPECT_EQ(ix,0);
    EXPECT_EQ(iy,0);

    // requesting a point right/above the domain should clamp to the last cell
    bilinear(x,y,3,3,50.0,50.0,ix,iy,coef);
    EXPECT_EQ(ix,1);
    EXPECT_EQ(iy,1);
}
