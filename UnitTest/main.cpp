//
//  main.cpp
//  UnitTest
//
//  Created by Vladislav Metel on 01/09/2018.
//  Copyright © 2018 Владислав Метель. All rights reserved.
//

#include <iostream>
#include <gtest/gtest.h>


TEST (SquareRootTest, ZeroAndNegativeNos) {
    ASSERT_EQ (0.0, 0.0);
    ASSERT_EQ (-1, -1);
}

int main(int argc, char * argv[]) {
    // insert code here...
    std::cout << "Hello, World!\n";
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
