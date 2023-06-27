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
    {
        LinearAlgebra::Mat<3, 3> in1({ { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } });
    }
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

consteval void testLinearAlgebra()
{
    static_assert(testVecOps(), "Failed Vector operations");
    static_assert(testPosOps(), "Failed Position operations");
    static_assert(testRayOps(), "Failed Ray operations");
}
