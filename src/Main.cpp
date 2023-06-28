#include "LinearAlgebra.hpp"

consteval void testLinearAlgebra();

int main()
{
    if consteval {
        testLinearAlgebra();
    }

    return 0;
}

consteval bool testVecOps()
{
    using namespace LinearAlgebra;
    bool has_passed = true;
    { // comparisons
        Vec<10> vec1 { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
        Vec<10> vec2 { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
        Vec<10> vec3 { 1, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
        Vec<9> vec4 { 0, 1, 2, 3, 4, 5, 6, 7, 8 };

        has_passed &= vec1 == vec2;
        has_passed &= vec2 != vec3;
        has_passed &= vec3 != vec4;
    }
    { // vector add, sub, mult, divide & indexing.
        Vec<10> in1 { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
        Vec<10> in2 { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };

        Vec<10> out1 { 1, 3, 5, 7, 9, 11, 13, 15, 17, 19 };
        Vec<10> out2 { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
        Vec<10> out3 { 0, 2, 6, 12, 20, 30, 42, 56, 72, 90 };
        Vec<10> out4 { 0, 0.5, 0.666667, 0.75, 0.8, 0.833333, 0.857143, 0.875, 0.888889, 0.9 };

        has_passed &= out1 == (in1 + in2);
        has_passed &= out2 == (in1 - in2);
        has_passed &= out3 == (in1 * in2);
        has_passed &= out4 == (in1 / in2);
        has_passed &= 7 == in1[7];
    }
    { // scalar mult, divide
        Vec<10> in1 { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };

        Vec<10> out1 { 2, 4, 6, 8, 10, 12, 14, 16, 18, 20 };
        Vec<10> out2 { 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5 };

        has_passed &= out1 == (in1 * 2);
        has_passed &= out1 == (2 * in1);
        has_passed &= out2 == (in1 / 2);
    }

    { // normalising, length and lengthSquared
        Vec<3> in1 { 3, 4, 0 };
        Vec<3> out1 { 0.6f, 0.8f, 0.0f };

        has_passed &= in1.getLength() == 5.f;
        has_passed &= in1.getLengthSquared() == 25.f;
        has_passed &= out1 == in1.getNormalised();
    }
    return has_passed;
}

consteval bool testPosOps()
{
    using namespace LinearAlgebra;
    bool has_passed = true;
    {
        Pos<3, float> in1 { 1, 2, 3 };
        Pos<3, float> in2 { 1, 1, 3 };

        Vec<3, float> out1 { 1, 2, 3 };
        has_passed &= out1 == static_cast<Vec<3, float>>(in1);
        has_passed &= in1 == static_cast<Pos<3, float>>(out1);
        has_passed &= 1.0f == LinearAlgebra::distance(in1, in2);
    }
    return has_passed;
}

consteval bool testMatOps()
{
    using namespace LinearAlgebra;
    bool has_passed = true;

    // Initializer list constructor and access
    Mat<3, 3> in1({ { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } });
    Mat<3, 3> in2({ { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } });
    Mat<3, 3> in3({ { 1, 2, 4 }, { 4, 5, 6 }, { 7, 8, 9 } });

    // Equality checks
    has_passed &= in1 == in2;
    has_passed &= in1 != in3;

    // Scalar multiplication
    Mat<3, 3> out1({ { 2, 4, 6 }, { 8, 10, 12 }, { 14, 16, 18 } });
    has_passed &= out1 == (in1 * 2);

    // Matrix-matrix multiplication
    Mat<3, 3> in4({ { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } });
    Mat<3, 3> out2({ { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } });
    has_passed &= out2 == (in1 * in4);

    // Identity matrix multiplication
    Mat<2, 2> id2({ { 1, 0 }, { 0, 1 } });
    Mat<2, 2> rand2({ { 3, 2 }, { 1, 4 } });
    has_passed &= rand2 == (id2 * rand2);
    has_passed &= rand2 == (rand2 * id2);

    // More complex multiplication
    Mat<3, 3> m1({ { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } });
    Mat<3, 3> m2({ { 10, 11, 12 }, { 13, 14, 15 }, { 16, 17, 18 } });
    Mat<3, 3> m3({ { 84, 90, 96 }, { 201, 216, 231 }, { 318, 342, 366 } });
    has_passed &= m3 == (m1 * m2);

    // Empty constructor
    Mat<3, 3> empty;
    has_passed &= empty[0][0] == 0 && empty[2][2] == 0; // checking some elements to be zero

    // Size inequality check
    Mat<2, 3> in5({ { 1, 2, 3 }, { 4, 5, 6 } });
    has_passed &= in5 != in1;

    return has_passed;
}

consteval bool testRayOps()
{
    using namespace LinearAlgebra;
    bool has_passed = true;
    {
        LinearAlgebra::Pos<3, float> origin = { 1, 2, 3 };
        LinearAlgebra::Vec<3, float> direction = { 1, 1, 1 };
        LinearAlgebra::Ray<3, float> ray = { origin, direction };
        LinearAlgebra::Pos<3, float> expected = { 2.154701076f, 3.154701076f, 4.154701076f };
        has_passed &= origin == ray.getOrigin();
        has_passed &= direction.getNormalised() == ray.getDirection();
        has_passed &= expected == ray.getPointAlongRay(2);
    }
    return has_passed;
}

consteval bool testMatRotOps()
{
    using namespace LinearAlgebra;
    bool has_passed = true;

    // Test for zero rotation
    {
        Mat<3, 3> rotMat = getRotationMat3x3<float>(0.0, 0.0, 0.0);
        Vec<3> vec(1, 0, 0);
        Vec<3> result = dotProduct(rotMat, vec);
        has_passed &= result == vec;
    }

    // Test for a 90-degree rotation around the Z-axis
    {
        Mat<3, 3> rotMat = getRotationMat3x3<float>(0.0, 0.0, toRadians(90));
        Vec<3> vec(1, 0, 0);
        Vec<3> expected_result(0, 1, 0);
        Vec<3> result = dotProduct(rotMat, vec);
        has_passed &= result == expected_result;
    }

    // Test for a 180-degree rotation around the Y-axis
    {
        Mat<3, 3, float> rotMat = getRotationMat3x3<float>(0.0, toRadians(180), 0.0);
        Vec<3, float> vec(1, 0, 0);
        Vec<3, float> expected_result(-1, 0, 0);
        Vec<3, float> result = dotProduct(rotMat, vec);
        has_passed &= result == expected_result;
    }

    return has_passed;
}

consteval bool testMatVecOps()
{
    namespace LA = LinearAlgebra;
    bool has_passed = true;
    {
        LA::Mat<3, 3, float> mat({ { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } });
        LA::Vec<3, float> vec { 1, 2, 3 };
        LA::Vec<3, float> mat_vec { 14, 32, 50 };
        LA::Vec<3, float> vec_mat { 30, 36, 42 };

        has_passed &= mat_vec == LA::dotProduct(mat, vec);
        has_passed &= vec_mat == LA::dotProduct(vec, mat);
    }

    return has_passed;
}

consteval void testLinearAlgebra()
{
    static_assert(testVecOps(), "Failed Vector operations");
    static_assert(testPosOps(), "Failed Position operations");
    static_assert(testRayOps(), "Failed Ray operations");
    static_assert(testMatOps(), "Failed Matrix operations");
    static_assert(testMatVecOps(), "Failed Matrix and Vector operations");
    static_assert(testMatRotOps(), "Failed Matrix rotation operations");
}
