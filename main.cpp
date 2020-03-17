#include <algorithm>
#include <array>
#include <stdarg.h>
#include <stdio.h>
#include <string>
#include <vector>

//-----------------------------------------

static const char* c_outFileName = "out.csv";
static const size_t c_numIterations = 25;

//-----------------------------------------

using Vec3 = std::array<float, 3>;

template <size_t N>
const std::array<float, N> operator + (const std::array<float, N>& A, const std::array<float, N>& B)
{
    std::array<float, N> ret;
    for (size_t i = 0; i < N; ++i)
        ret[i] = A[i] + B[i];
    return ret;
}

template <size_t N>
const std::array<float, N> operator - (const std::array<float, N>& A, const std::array<float, N>& B)
{
    std::array<float, N> ret;
    for (size_t i = 0; i < N; ++i)
        ret[i] = A[i] - B[i];
    return ret;
}

template <size_t N>
const std::array<float, N> operator * (const std::array<float, N>& A, float B)
{
    std::array<float, N> ret;
    for (size_t i = 0; i < N; ++i)
        ret[i] = A[i] * B;
    return ret;
}

//-----------------------------------------

struct CSVFile
{
    void SetCell(size_t x, size_t y, const char* format, ...)
    {
        char buffer[4096];
        va_list args;
        va_start(args, format);
        vsprintf_s(buffer, format, args);
        va_end(args);

        if (rows.size() <= y)
            rows.resize(y + 1);

        row& r = rows[y];
        if (r.size() <= x)
            r.resize(x + 1);

        r[x] = buffer;
    }

    void Save(const char* fileName)
    {
        FILE* file = nullptr;
        fopen_s(&file, fileName, "w+t");

        for (const row& r : rows)
        {
            bool first = true;
            for (const std::string& s : r)
            {
                fprintf(file, "%s\"%s\"", first ? "" : ",", s.c_str());
                first = false;
            }
            fprintf(file, "\n");
        }

        fclose(file);
    }

    using row = std::vector<std::string>;
    std::vector<row> rows;
};

struct TestDesc
{
    struct NewtonAndHalley
    {
        float initialGuess;
    };

    struct Secant
    {
        float initialGuess;
        float priorInitialGuess;
    };

    struct Bisection
    {
        float min;
        float max;
    };

    std::vector<NewtonAndHalley> newton;
    std::vector<NewtonAndHalley> halley;
    std::vector<Secant> secant;
    std::vector<Bisection> bisection;
};

template <typename LAMBDA_VALUE, typename LAMBDA_FIRST_DERIVATIVE>
void RootFind_Newton(float initialGuess, CSVFile& csv, const LAMBDA_VALUE& lambdaValue, const LAMBDA_FIRST_DERIVATIVE& lambdaFirstDerivative)
{
    size_t columnIndex = csv.rows[0].size();

    csv.SetCell(columnIndex + 0, 0, "Newton X  [x=%0.2f]", initialGuess);
    csv.SetCell(columnIndex + 1, 0, "Newton Y  [x=%0.2f]", initialGuess);
    csv.SetCell(columnIndex + 2, 0, "Newton Abs Y  [x=%0.2f]", initialGuess);

    float x = initialGuess;
    for (size_t iterationIndex = 1; iterationIndex <= c_numIterations; ++iterationIndex)
    {
        float y = lambdaValue(x);
        float yPrime = lambdaFirstDerivative(x);

        csv.SetCell(columnIndex + 0, iterationIndex, "%f", x);
        csv.SetCell(columnIndex + 1, iterationIndex, "%f", y);
        csv.SetCell(columnIndex + 2, iterationIndex, "%f", abs(y));

        // If we will divide by zero when updating x, we are as close as we are going to get.
        // We should return out of the function, but just doing a continue so the csv isn't missing data.
        if (yPrime == 0.0f)
            continue;

        x = x - y / yPrime;
    }
}

template <typename LAMBDA_VALUE, typename LAMBDA_FIRST_DERIVATIVE, typename LAMBDA_SECOND_DERIVATIVE>
void RootFind_Halley(
    float initialGuess,
    CSVFile& csv,
    const LAMBDA_VALUE& lambdaValue,
    const LAMBDA_FIRST_DERIVATIVE& lambdaFirstDerivative,
    const LAMBDA_SECOND_DERIVATIVE& lambdaSecondDerivative
)
{
    size_t columnIndex = csv.rows[0].size();

    csv.SetCell(columnIndex + 0, 0, "Halley X  [x=%0.2f] ", initialGuess);
    csv.SetCell(columnIndex + 1, 0, "Halley Y  [x=%0.2f]", initialGuess);
    csv.SetCell(columnIndex + 2, 0, "Halley Abs Y  [x=%0.2f]", initialGuess);

    float x = initialGuess;
    for (size_t iterationIndex = 1; iterationIndex <= c_numIterations; ++iterationIndex)
    {
        float y = lambdaValue(x);
        float yPrime = lambdaFirstDerivative(x);
        float yPrimePrime = lambdaSecondDerivative(x);

        csv.SetCell(columnIndex + 0, iterationIndex, "%f", x);
        csv.SetCell(columnIndex + 1, iterationIndex, "%f", y);
        csv.SetCell(columnIndex + 2, iterationIndex, "%f", abs(y));

        x = x - (2.0f * y * yPrime) / (2.0f * yPrime * yPrime - y * yPrimePrime);
    }
}

template <typename LAMBDA_VALUE>
void RootFind_Secant(float initialGuess, float priorInitialGuess, CSVFile& csv, const LAMBDA_VALUE& lambdaValue)
{
    size_t columnIndex = csv.rows[0].size();

    csv.SetCell(columnIndex + 0, 0, "Secant X  [x=%0.2f, px=%0.2f]", initialGuess, priorInitialGuess);
    csv.SetCell(columnIndex + 1, 0, "Secant Y  [x=%0.2f, px=%0.2f]", initialGuess, priorInitialGuess);
    csv.SetCell(columnIndex + 2, 0, "Secant Abs Y  [x=%0.2f, px=%0.2f]", initialGuess, priorInitialGuess);

    float priorX = priorInitialGuess;
    float priorY = lambdaValue(priorX);
    float x = initialGuess;
    for (size_t iterationIndex = 1; iterationIndex <= c_numIterations; ++iterationIndex)
    {
        float y = lambdaValue(x);
        float yPrime = (y - priorY) / (x - priorX);  // Secant is just newton using finite differences for 1st derivative!

        csv.SetCell(columnIndex + 0, iterationIndex, "%f", x);
        csv.SetCell(columnIndex + 1, iterationIndex, "%f", y);
        csv.SetCell(columnIndex + 2, iterationIndex, "%f", abs(y));

        // If yPrime divided by zero, or we will divide by zero when updating x, we are as close as we are going to get.
        // We should return out of the function, but just doing a continue so the csv isn't missing data.
        if ((x - priorX) == 0.0f || yPrime == 0.0f)
            continue;

        priorX = x;
        priorY = y;
        x = x - y / yPrime;
    }
}

template <typename LAMBDA_VALUE>
void RootFind_Bisection(float minX, float maxX, CSVFile& csv, const LAMBDA_VALUE& lambdaValue)
{
    size_t columnIndex = csv.rows[0].size();

    csv.SetCell(columnIndex + 0, 0, "Bisect X  [%0.2f, %0.2f]", minX, maxX);
    csv.SetCell(columnIndex + 1, 0, "Bisect Y  [%0.2f, %0.2f]", minX, maxX);
    csv.SetCell(columnIndex + 2, 0, "Bisect Abs Y  [%0.2f, %0.2f]", minX, maxX);

    float minY = lambdaValue(minX);
    float maxY = lambdaValue(maxX);

    if (minY > maxY)
    {
        std::swap(minX, maxX);
        std::swap(minY, maxY);
    }

    // y signs need to be opposite
    if (minY > 0.0f || maxY < 0.0f)
        return;

    for (size_t iterationIndex = 1; iterationIndex <= c_numIterations; ++iterationIndex)
    {
        float midX = (minX + maxX) / 2.0f;
        float midY = lambdaValue(midX);

        csv.SetCell(columnIndex + 0, iterationIndex, "%f", midX);
        csv.SetCell(columnIndex + 1, iterationIndex, "%f", midY);
        csv.SetCell(columnIndex + 2, iterationIndex, "%f", abs(midY));

        if (midY < 0.0f)
        {
            minX = midX;
            minY = midY;
        }
        else
        {
            maxX = midX;
            maxY = midY;
        }
    }
}

template <typename LAMBDA_VALUE, typename LAMBDA_FIRST_DERIVATIVE, typename LAMBDA_SECOND_DERIVATIVE>
void DoTests(
    const TestDesc& descs,
    CSVFile& csv,
    const LAMBDA_VALUE& lambdaValue,
    const LAMBDA_FIRST_DERIVATIVE& lambdaFirstDerivative, 
    const LAMBDA_SECOND_DERIVATIVE& lambdaSecondDerivative
)
{
    for (const auto& desc : descs.newton)
        RootFind_Newton(desc.initialGuess, csv, lambdaValue, lambdaFirstDerivative);

    for (const auto& desc : descs.halley)
        RootFind_Halley(desc.initialGuess, csv, lambdaValue, lambdaFirstDerivative, lambdaSecondDerivative);

    for (const auto& desc : descs.secant)
        RootFind_Secant(desc.initialGuess, desc.priorInitialGuess, csv, lambdaValue);

    for (const auto& desc : descs.bisection)
        RootFind_Bisection(desc.min, desc.max, csv, lambdaValue);
}

void RayVsSphereTest(CSVFile& csv)
{
    // put a blank column, then a label column before the tests
    csv.SetCell(csv.rows[0].size(), 0, "");
    csv.SetCell(csv.rows[0].size(), 0, "Ray Vs Sphere");

    static const Vec3 c_rayPos = { 0.0f, 0.0f, 0.0f };
    static const Vec3 c_rayDir = { 0.0f, 0.0f, 1.0f };

    static const Vec3 c_spherePos = { 0.2f, 0.3f, 5.0f };
    static const float c_sphereRadius = 2.0f;

    static const Vec3 c_relativePos = c_rayPos - c_spherePos;  // make the sphere center be the origin

    // For a value of t, return the squared distance of the ray to the sphere.
    // Squared distance because derivatives are easier without the square root.
    // SquaredDistance = Magnitude(RelativePos + RayDir * T) - SphereRadius^2
    auto Value = [] (float t) -> float
    {
        Vec3 diff = c_relativePos + c_rayDir * t;
        return diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2] - c_sphereRadius * c_sphereRadius;
    };

    // For a value of t, return the first derivative of the squared distance of the ray to the sphere.
    // SquaredDistance' = 2 * RelativePos.x * RayDir.x + RayDir.x^2 * T 
    // Add in the .y and .z versions too.
    auto FirstDerivative = [&](float t) -> float
    {
#if 1
        float ret = 0.0f;
        for (size_t index = 0; index < 3; ++index)
            ret += 2.0f * c_relativePos[index] * c_rayDir[index] + 2.0f * c_rayDir[index] * c_rayDir[index] * t;
        return ret;
#else
        // central difference method for calculating derivative
        static const float c_epsilon = 0.01f;
        float numericalDerivative = (Value(t + c_epsilon) - Value(t - c_epsilon)) / (2.0f * c_epsilon);
        return numericalDerivative;
#endif
    };

    // For a value of t, return the second derivative of the squared distance of the ray to the sphere.
    // SquaredDistance'' = 2 * RayDir.x^2
    // Add in the .y and .z versions too.
    auto SecondDerivative = [&](float t)
    {
#if 1
        float ret = 0.0f;
        for (size_t index = 0; index < 3; ++index)
            ret += 2.0f * c_rayDir[index] * c_rayDir[index];
        return ret;
#else
        // central difference method for calculating derivative
        static const float c_epsilon = 0.01f;
        float numericalDerivative = (FirstDerivative(t + c_epsilon) - FirstDerivative(t - c_epsilon)) / (2.0f * c_epsilon);
        return numericalDerivative;
#endif
    };

    TestDesc desc;
    desc.newton = { {0.0f} };
    desc.halley = desc.newton;
    desc.secant = { {0.1f, 0.2f} };
    desc.bisection = { { 0.0f, 6.0f} };
    DoTests(desc, csv,
        Value,
        FirstDerivative,
        SecondDerivative
    );
}

int main(int argc, char** argv)
{
    CSVFile csv;

    // put iteration count down the left column
    csv.SetCell(0, 0, "");
    for (size_t index = 1; index <= c_numIterations; ++index)
        csv.SetCell(0, index, "%zu", index);

    // y = x^2-1
    // y' = 2x
    // y'' = 2
    {
        // put a blank column, then a label column before the tests
        csv.SetCell(csv.rows[0].size(), 0, "");
        csv.SetCell(csv.rows[0].size(), 0, "x^2-1");

        TestDesc desc;
        desc.newton = { {0.5f}, {10.0f} };
        desc.halley = desc.newton;
        desc.secant = { {0.5f, 0.4f}, {10.0f, 9.0f} };
        desc.bisection = { {-1.5f, 0.5f}, {-100.0f, 0.5f} };
        DoTests(desc, csv,
            [](float x) { return x * x - 1.0f; },
            [](float x) { return 2.0f * x; },
            [](float x) { return 2.0f; }
        );
    }

    // y = sin(x)
    // y' = cos(x)
    // y'' = -sin(x)
    {
        // put a blank column, then a label column before the tests
        csv.SetCell(csv.rows[0].size(), 0, "");
        csv.SetCell(csv.rows[0].size(), 0, "sin(x)");

        TestDesc desc;
        desc.newton = { {6.0f} };
        desc.halley = desc.newton;
        desc.secant = { {6.0f, 5.9f} };
        desc.bisection = { {3.5f, 7.0f} };
        DoTests(desc, csv,
            [](float x) { return sinf(x); },
            [](float x) { return cosf(x); },
            [](float x) { return -sinf(x); }
        );
    }

    // y = x^2-x-1
    // y' = 2x-1
    // y'' = 2
    // AKA calculate the golden ratio
    {
        // put a blank column, then a label column before the tests
        csv.SetCell(csv.rows[0].size(), 0, "");
        csv.SetCell(csv.rows[0].size(), 0, "x^2-x-1");

        TestDesc desc;
        desc.newton = { {0.5f} };
        desc.halley = desc.newton;
        desc.secant = { {0.5f, 0.4f} };
        desc.bisection = { {-1.5f, 0.5f} };
        DoTests(desc, csv,
            [](float x) { return x * x - x - 1.0f; },
            [](float x) { return 2.0f * x - 1.0f; },
            [](float x) { return 2.0f; }
        );
    }

    // y = x^4-3x^3-20x-100
    // y' = 4x^3-9x^2-20
    // y'' = 12x^2-18x
    // AKA calculate the golden ratio
    {
        // put a blank column, then a label column before the tests
        csv.SetCell(csv.rows[0].size(), 0, "");
        csv.SetCell(csv.rows[0].size(), 0, "x^4-3x^3-20x-100");

        TestDesc desc;
        desc.newton = { {0.5f} };
        desc.halley = desc.newton;
        desc.secant = { {0.5f, 0.4f} };
        desc.bisection = { {-2.5f, -1.5f} };
        DoTests(desc, csv,
            [](float x) { return (x*x*x*x)-3.0f*(x*x*x) - 20.0f*x - 100.0f; },
            [](float x) { return 4.0f * (x * x * x) - 9.0f * (x * x) - 20.0f; },
            [](float x) { return 12.0f * (x * x) - 18.0f * x; }
        );
    }

    // ERROR CASE: This only has imaginary roots. Also bisect bounds are not met
    // y = x^2+1
    // y' = 2x
    // y'' = 2
    {
        // put a blank column, then a label column before the tests
        csv.SetCell(csv.rows[0].size(), 0, "");
        csv.SetCell(csv.rows[0].size(), 0, "x^2+1");

        TestDesc desc;
        desc.newton = { {0.5f} };
        desc.halley = desc.newton;
        desc.secant = { {0.5f, 0.4f} };
        desc.bisection = { {-1.5f, 0.5f} };
        DoTests(desc, csv,
            [](float x) { return x * x + 1.0f; },
            [](float x) { return 2.0f * x; },
            [](float x) { return 2.0f; }
        );
    }

    // Ray vs Sphere root finding
    RayVsSphereTest(csv);

    // save the results
    csv.Save(c_outFileName);
    return 0;
}

/*

NOTES:
* look at your email but also...
 * secant needs: f(x), an initial guess, and a prior initial guess
  * secant is just newton, but using finite differences to get f'(x)
 * newton needs: f(x), f'(x) and an initial guess.
 * halley needs: f(x), f'(x), f''(x) and an initial guess
 * bisection needs: f(x), min and max, bounding the zero, and having f(min) and f(max) having opposite signs
  * could just go downhill if they didn't have opposite signs but that'd just be good for finding local minima (optimizing), not zeroes.

* show when there's only imaginary roots

* multi dimensional root finding - you can extend newton and the others easily enough. Halley getting 2nd derivative would involve a jacobian matrix maybe?
 * ray vs sphere is a simple example. can compare vs numerical method. ? does it take a single newton step since it's linear? so don't even need numerical solution apparently.
 * with a zero second derivative, halley decays to newton. so, lying and saying "0" is the same as just doing newton. Same if it really is 0.
 * ray vs sphere could be ray vs other mathematical shape. Problem is in coming up w/ derivatives. Can use Secant method with numerical derivatives though, or dual numbers.
 * we are doing numerical derivatives, which means newton decays to secant etc.
 * This is sort what we do in sphere marching: https://www.iquilezles.org/www/articles/distance/distance.htm
  * we get the length of the gradient, and the value at that point, and know that nothing can be closer than how long it would take for the value to reach zero along that gradient.

- if you know the root to a polynomial, can you simplify it?

- if you know a root to a polynomial, can you simplify it for the remaining roots?


! put derivation for sphere thing on blog post.

! amazing how well newton does for just a few steps.

*/