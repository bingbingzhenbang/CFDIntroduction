#include "SubsonicSupersonicIsentropicNozzleFlow.h"
#include "PurelySubsonicIsentropicNozzleFlow.h"
#include <fstream>

int main()
{
    std::ofstream fout1("SubsonicSupersonicIsentropicNozzleFlow.csv");
    SubsonicSupersonic solve1(0.1, 30);
    solve1.Calculate(1.4, 0.5, 1400, fout1);
	return 0;
}