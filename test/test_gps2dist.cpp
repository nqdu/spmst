#include "shared/gps2dist.hpp"
#include <gtest/gtest.h>
#include <cmath>

// note the unusual argument order: gps2dist(lon0,lon1,lat0,lat1,radius)

TEST(Gps2dist, SamePointIsZero)
{
    EXPECT_NEAR(gps2dist(10.0,10.0,20.0,20.0,6371.0), 0.0, 1e-3);
}

TEST(Gps2dist, QuarterOfEquator)
{
    // (0,0) to (90,0) along the equator is a quarter of the great circle
    float expect = M_PI / 2.0 * 6371.0;
    EXPECT_NEAR(gps2dist(0.0,90.0,0.0,0.0,6371.0), expect, 1.0);
}

TEST(Gps2dist, PoleToPole)
{
    // north pole to south pole is half the great circle, regardless of longitude
    float expect = M_PI * 6371.0;
    EXPECT_NEAR(gps2dist(0.0,0.0,90.0,-90.0,6371.0), expect, 1.0);
}

TEST(Gps2dist, OneDegreeLatitudeNearEquator)
{
    // one degree of latitude is about 111.19 km on Earth's surface
    EXPECT_NEAR(gps2dist(0.0,0.0,0.0,1.0,6371.0), 111.19, 0.01);
}

TEST(Gps2dist, Symmetric)
{
    float d1 = gps2dist(10.0,20.0,30.0,40.0,6371.0);
    float d2 = gps2dist(20.0,10.0,40.0,30.0,6371.0);
    EXPECT_NEAR(d1, d2, 1e-3);
}
