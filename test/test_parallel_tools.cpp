#include "shared/parallel_tools.hpp"
#include <gtest/gtest.h>
#include <vector>

TEST(AllocateTasks, EvenSplit)
{
    // 10 tasks, 2 procs -> 5 each
    int start,end;
    allocate_tasks(10,2,0,start,end);
    EXPECT_EQ(start,0); EXPECT_EQ(end,4);
    allocate_tasks(10,2,1,start,end);
    EXPECT_EQ(start,5); EXPECT_EQ(end,9);
}

TEST(AllocateTasks, UnevenSplitGivesExtraToFirstRanks)
{
    // 10 tasks, 3 procs -> sizes 4,3,3
    int start,end;
    allocate_tasks(10,3,0,start,end);
    EXPECT_EQ(start,0); EXPECT_EQ(end,3);
    allocate_tasks(10,3,1,start,end);
    EXPECT_EQ(start,4); EXPECT_EQ(end,6);
    allocate_tasks(10,3,2,start,end);
    EXPECT_EQ(start,7); EXPECT_EQ(end,9);
}

TEST(AllocateTasks, AllTasksCoveredExactlyOnce)
{
    int ntasks = 17, nprocs = 5;
    std::vector<int> count(ntasks,0);
    for(int rank = 0; rank < nprocs; rank++){
        int start,end;
        allocate_tasks(ntasks,nprocs,rank,start,end);
        for(int i = start; i <= end; i++) count[i] += 1;
    }
    for(int i = 0; i < ntasks; i++){
        EXPECT_EQ(count[i],1) << "task " << i << " covered " << count[i] << " times";
    }
}

TEST(AllocateTasks, MoreProcsThanTasksGivesEmptyRange)
{
    // more procs than tasks -> ranks beyond ntasks get an empty (start>end) range
    int start,end;
    allocate_tasks(3,5,4,start,end);
    EXPECT_GT(start,end);
}
