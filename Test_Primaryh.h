#pragma once
#include <stdexcept>
#include "Primary.h"

using namespace std;
bool testPrimary(Primary& distribution)
{
    double* array = new double[4];
    int array_of_nu[]{ 1, 2, 3, 4, 5, 6, 7 };
    double array_of_disp[]{ 2, 4, 6, 8, 10, 12, 14 };
    double array_of_exc[]{ 3, 1.5, 1, 0.75, 0.6, 0.5, 0.429 };
    double array_of_densities[]{ 0.5 , 0.25, 0.188, 0.156, 0.137, 0.123, 0.113 };
    for (int i = 0; i < 7; i++)
        try
    {
        distribution.set_n(array_of_nu[i]);
        array = distribution.get_metrics();
        if (abs(array[1] - array_of_disp[i]) > 0.01)
        {
            throw 0;
        }
        if (abs(array[3] - array_of_exc[i]) > 0.01)
        {
            throw 0;
        }
        if (abs(distribution.Laplace_function(0) - array_of_densities[i]) > 0.01)
        {
            throw 0;
        }
        
    }
    catch (int error)
    {
        if (error == 0)
        {
            cout << "Ошибка при n =" << distribution.get_mu() << endl;
            return false;
        }
    }
    return true;
}