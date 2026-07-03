#include "invparam.hpp"
#include <gtest/gtest.h>
#include <fstream>
#include <cstdio>

class InvParamTest : public ::testing::Test {
protected:
    std::string filename = "test_invparam_tmpfile.in";

    void write_file(const std::string &content)
    {
        std::ofstream out(filename);
        out << content;
        out.close();
    }

    void TearDown() override
    {
        std::remove(filename.c_str());
    }
};

TEST_F(InvParamTest, ParsesAllFields)
{
    write_file(
        "# SPMST inversion Parameters\n"
        "NITERS = 6\n"
        "ITER_CURRENT = 2\n"
        "MIN_VELOC = 2.0\n"
        "MAX_VELOC = 5.0\n"
        "SPHERICAL = 1\n"
        "SYN_TEST = 1\n"
        "NOISE_LEVEL = 0.1\n"
        "SMOOTH = 2.0\n"
        "DAMP = 0.01\n"
        "NTHREADS = 2\n"
    );

    InverseParamsBase param;
    param.read_file(filename.c_str());

    EXPECT_EQ(param.maxiter,6);
    EXPECT_EQ(param.iter_cur,2);
    EXPECT_FLOAT_EQ(param.minvel,2.0f);
    EXPECT_FLOAT_EQ(param.maxvel,5.0f);
    EXPECT_TRUE(param.is_spherical);
    EXPECT_EQ(param.ifsyn,1);
    EXPECT_FLOAT_EQ(param.noiselevel,0.1f);
    EXPECT_FLOAT_EQ(param.smooth,2.0f);
    EXPECT_FLOAT_EQ(param.damp,0.01f);
    EXPECT_EQ(param.nthreads,2);
}

TEST_F(InvParamTest, DefaultsIterCurrentToZeroWhenMissing)
{
    write_file(
        "NITERS = 6\n"
        "MIN_VELOC = 2.0\n"
        "MAX_VELOC = 5.0\n"
        "SPHERICAL = 0\n"
        "SYN_TEST = 0\n"
        "SMOOTH = 2.0\n"
        "DAMP = 0.01\n"
        "NTHREADS = 2\n"
    );

    InverseParamsBase param;
    param.read_file(filename.c_str());

    EXPECT_EQ(param.iter_cur,0);
    EXPECT_FALSE(param.is_spherical);
}
