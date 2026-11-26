#include <iostream>

#include "FiniteTempRelaxer.h"

int main()
{
    FiniteTempRelaxer relaxer;
    relaxer.read_parameters();

    const bool ok = relaxer.run();
    if (!ok)
    {
        std::cout << "Relaxation did not converge.\n";
        return 1;
    }

    std::cout << "Relaxation complete.\n";
    return 0;
}