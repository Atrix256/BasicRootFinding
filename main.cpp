#include <algorithm>
#include <stdarg.h>
#include <stdio.h>
#include <string>
#include <vector>

static const char* c_outFileName = "out.csv";
static const size_t c_numIterations = 25;

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
    const char* label = nullptr;

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

    std::vector<NewtonAndHalley> newtonAndHalley;
    std::vector<Secant> secant;
    std::vector<Bisection> bisection;
};

template <typename LAMBDA_VALUE, typename LAMBDA_FIRST_DERIVATIVE>
void RootFind_Newton(const char* label, float initialGuess, CSVFile& csv, const LAMBDA_VALUE& lambdaValue, const LAMBDA_FIRST_DERIVATIVE& lambdaFirstDerivative)
{
    size_t columnIndex = csv.rows[0].size();

    csv.SetCell(columnIndex + 0, 0, "Newton X  %s [x=%0.2f]", label, initialGuess);
    csv.SetCell(columnIndex + 1, 0, "Newton Y  %s [x=%0.2f]", label, initialGuess);

    float x = initialGuess;
    for (size_t iterationIndex = 1; iterationIndex <= c_numIterations; ++iterationIndex)
    {
        float y = lambdaValue(x);
        float yPrime = lambdaFirstDerivative(x);

        csv.SetCell(columnIndex + 0, iterationIndex, "%f", x);
        csv.SetCell(columnIndex + 1, iterationIndex, "%f", y);

        x = x - y / yPrime;
    }
}

template <typename LAMBDA_VALUE, typename LAMBDA_FIRST_DERIVATIVE, typename LAMBDA_SECOND_DERIVATIVE>
void RootFind_Halley(
    const char* label,
    float initialGuess,
    CSVFile& csv,
    const LAMBDA_VALUE& lambdaValue,
    const LAMBDA_FIRST_DERIVATIVE& lambdaFirstDerivative,
    const LAMBDA_SECOND_DERIVATIVE& lambdaSecondDerivative
)
{
    size_t columnIndex = csv.rows[0].size();

    csv.SetCell(columnIndex + 0, 0, "Halley X  %s [x=%0.2f] ", label, initialGuess);
    csv.SetCell(columnIndex + 1, 0, "Halley Y  %s [x=%0.2f]", label, initialGuess);

    float x = initialGuess;
    for (size_t iterationIndex = 1; iterationIndex <= c_numIterations; ++iterationIndex)
    {
        float y = lambdaValue(x);
        float yPrime = lambdaFirstDerivative(x);
        float yPrimePrime = lambdaSecondDerivative(x);

        csv.SetCell(columnIndex + 0, iterationIndex, "%f", x);
        csv.SetCell(columnIndex + 1, iterationIndex, "%f", y);

        x = x - (2.0f * y * yPrime) / (2.0f * yPrime * yPrime - y * yPrimePrime);
    }
}

template <typename LAMBDA_VALUE>
void RootFind_Secant(const char* label, float initialGuess, float priorInitialGuess, CSVFile& csv, const LAMBDA_VALUE& lambdaValue)
{
    size_t columnIndex = csv.rows[0].size();

    csv.SetCell(columnIndex + 0, 0, "Secant X  %s [x=%0.2f, px=%0.2f]", label, initialGuess, priorInitialGuess);
    csv.SetCell(columnIndex + 1, 0, "Secant Y  %s [x=%0.2f, px=%0.2f]", label, initialGuess, priorInitialGuess);

    float priorX = priorInitialGuess;
    float priorY = lambdaValue(priorX);
    float x = initialGuess;
    for (size_t iterationIndex = 1; iterationIndex <= c_numIterations; ++iterationIndex)
    {
        float y = lambdaValue(x);
        float yPrime = (y - priorY) / (x - priorX);  // Secant is just newton using finite differences for 1st derivative!

        csv.SetCell(columnIndex + 0, iterationIndex, "%f", x);
        csv.SetCell(columnIndex + 1, iterationIndex, "%f", y);

        priorX = x;
        priorY = y;
        x = x - y / yPrime;
    }
}

template <typename LAMBDA_VALUE>
void RootFind_Bisection(const char* label, float minX, float maxX, CSVFile& csv, const LAMBDA_VALUE& lambdaValue)
{
    size_t columnIndex = csv.rows[0].size();

    csv.SetCell(columnIndex + 0, 0, "Bisect X  %s [%0.2f, %0.2f]", label, minX);
    csv.SetCell(columnIndex + 1, 0, "Bisect Y  %s [%0.2f, %0.2f]", label, maxX);

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
        float midX = (min + max) / 2.0f;
        float midY = lambdaValue(mid);

        float yPrime = (y - priorY) / (x - priorX);  // Secant is just newton using finite differences for 1st derivative!

        csv.SetCell(columnIndex + 0, iterationIndex, "%f", x);
        csv.SetCell(columnIndex + 1, iterationIndex, "%f", y);

        x = x - y / yPrime;
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
    for (const auto& desc : descs.newtonAndHalley)
    {
        RootFind_Newton(descs.label, desc.initialGuess, csv, lambdaValue, lambdaFirstDerivative);
        RootFind_Halley(descs.label, desc.initialGuess, csv, lambdaValue, lambdaFirstDerivative, lambdaSecondDerivative);
    }

    for (const auto& desc : descs.secant)
    {
        RootFind_Secant(descs.label, desc.initialGuess, desc.priorInitialGuess, csv, lambdaValue);
    }

    for (const auto& desc : descs.bisection)
    {
        RootFind_Bisection(descs.label, desc.min, desc.max, csv, lambdaValue);
    }
}

int main(int argc, char** argv)
{
    CSVFile csv;

    // put iteration count down the left column
    csv.SetCell(0, 0, "");
    for (size_t index = 1; index <= c_numIterations; ++index)
        csv.SetCell(0, index, "%zu", index);

    // y = x^2-1
    {
        TestDesc desc;
        desc.label = "y=x^2-1";
        desc.newtonAndHalley = { {0.5f}, {10.0f} };
        desc.secant = { {0.5f, 0.4f}, {10.0f, 9.0f} };
        desc.bisection = { {-1.5f, 0.5f}, {-100.0f, 0.5f} };
        DoTests(desc, csv,
            [](float x) { return x * x - 1.0f; },
            [](float x) { return 2.0f * x; },
            [](float x) { return 2.0f; }
        );
    }

    // put a space between this and the next test
    csv.SetCell(csv.rows[0].size(), 0, "");

    // save the results
    csv.Save(c_outFileName);
    return 0;
}

/*
TODO:

- different initial guesses?
- if you know the root to a polynomial, can you simplify it?
- for optimization, just do the first derivative of one of the other functions

* more functions to test. like could do ray vs sphere or something...


- could show what happens when there is no root (y=x^2+1 has imaginary roots)
- and what happens with bisection when it's needs are not met (opposite signs, bound the root)

NOTES:
* look at your email but also...
 * secant needs: f(x), an initial guess, and a prior initial guess
  * secant is just newton, but using finite differences to get f'(x)
 * newton needs: f(x), f'(x) and an initial guess.
 * halley needs: f(x), f'(x), f''(x) and an initial guess
 * bisection needs: f(x), min and max, bounding the zero, and having f(min) and f(max) having opposite signs
  * could just go downhill if they didn't have opposite signs but that'd just be good for finding local minima (optimizing), not zeroes.

*/