#pragma once

class PurelySubsonic
{
public:
    double GetA(double x)
    {
        double tmp = x - 1.5;
        if (0 <= x && x <= 1.5)
            return 1 + 2.2 * tmp * tmp;
        else if (x > 1.5 && x <= 3.0)
            return 1 + 0.2223 * tmp * tmp;
        else
            return 0.0;
    }
};