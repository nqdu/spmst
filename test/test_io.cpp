#include "shared/IO.hpp"
#include <gtest/gtest.h>
#include <fstream>
#include <cstdio>

class IOFileTest : public ::testing::Test {
protected:
    std::string filename = "test_io_tmpfile.txt";

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

TEST_F(IOFileTest, ReadFileParamSkipsCommentsAndBlankLines)
{
    write_file(
        "# this is a comment\n"
        "\n"
        "1 2.5\n"
        "hello\n"
    );

    std::ifstream fp(filename);
    int a; float b; std::string s;
    read_file_param(fp,a,b);
    read_file_param(fp,s);

    EXPECT_EQ(a,1);
    EXPECT_FLOAT_EQ(b,2.5f);
    EXPECT_EQ(s,"hello");
}

TEST_F(IOFileTest, ReadParRegexFindsNamedParam)
{
    write_file(
        "# header comment\n"
        "NITERS = 6\n"
        "MIN_VELOC = 2.0\n"
        "MAX_VELOC = 5.0\n"
    );

    std::ifstream fp(filename);
    int niters;
    float minvel,maxvel;
    EXPECT_EQ(read_par_regex("NITERS",niters,fp),0);
    EXPECT_EQ(niters,6);
    EXPECT_EQ(read_par_regex("MIN_VELOC",minvel,fp),0);
    EXPECT_FLOAT_EQ(minvel,2.0f);
    EXPECT_EQ(read_par_regex("MAX_VELOC",maxvel,fp),0);
    EXPECT_FLOAT_EQ(maxvel,5.0f);
}

TEST_F(IOFileTest, ReadParRegexMissingParamReturnsError)
{
    write_file("NITERS = 6\n");
    std::ifstream fp(filename);
    float damp;
    EXPECT_EQ(read_par_regex("DAMP",damp,fp),1);
}

TEST_F(IOFileTest, ReadParRegexIsOrderIndependent)
{
    write_file(
        "MAX_VELOC = 5.0\n"
        "NITERS = 6\n"
        "MIN_VELOC = 2.0\n"
    );

    std::ifstream fp(filename);
    int niters;
    // query MIN_VELOC first even though it appears last in the file
    float minvel;
    EXPECT_EQ(read_par_regex("MIN_VELOC",minvel,fp),0);
    EXPECT_FLOAT_EQ(minvel,2.0f);
    EXPECT_EQ(read_par_regex("NITERS",niters,fp),0);
    EXPECT_EQ(niters,6);
}
