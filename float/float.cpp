#include <cmath>
#include "BigFloat.h"
#include <iostream>
#include <bitset>

int main()
{
    /*double kak1 = 9;
    double kak2 = 3;

    int n = 1000000;
    BigFloat<3> s = 1;
    bool minus = true;
    BigFloat<3> j = 3;
    for (size_t i = 0; i < n; i++)
    {
        if (minus)
            s = s - (1 / j);
        else
            s = s + (1 / j);

        j = j + 2;
        minus = !minus;
    }
    s = s * 4;*/
    BigFloat<6> ts = 117827932;

    /* bigfloat<6> kaka = 10;
    ts = ts / kaka;*/

    std::string str = "1867.23984675230076000089042307840234328976666666666668796875611111";
    str = "2.0";
    BigFloat<5191> pi;
    try 
    {
        pi = BigFloat<5191>(str);
    }
    catch(const std::invalid_argument& e)
    {
        std::cerr << e.what() << std::endl;
    }
    catch (...)
    {
        std::cerr << "grosser error..?" << std::endl;
    }
    BigFloat<5191>::sqrt(pi);
    /*BigFloat<6> one(1);
    BigFloat<6> two(2);
    for (int i = 0; i < 129; i++)
    {
        pi = pi / two;
    }
    pi = pi + one;*/

    /*BigFloat<2> nob1 = kak1;      
    BigFloat<2> nob2 = kak2;
    BigFloat<2> c = nob1 / nob2;*/
    //

    std::vector<std::uint64_t> mant1 = { {0x9000000000000000}, {43} };
    int32_t exp1 = 3;
    bool sign = false;

    BigFloat<2> a(mant1, true, exp1);


    /*std::vector<std::uint64_t> mant2 = { {0x8000000000000000}, {0x0000000000000000} };
    int32_t exp2 = 3;

    BigFloat<2> b(mant2, true, exp2);

    BigFloat<2> c = nob1/nob2;


    for (auto i = 0; i <2; i++)
        std::cout << std::bitset<64> (mant1[i]);
    std::cout << std::endl;

    for (auto i = 0; i <2; i++)
        std::cout << std::bitset<64>(mant2[i]);
    std::cout << std::endl;*/
    //std::cout << std::bitset<64>(pi.getExp() + 127) << " dezimal: " << pi.getExp() << std::endl;
    //std::cout << std::bitset<11>(pi.getExp() + 1023);
    //std::cout << std::bitset<52>(pi.mantisse[pi.mantisse.size() - 1] >>(64-24));
    ////for ( size_t i = 1; i < pi.mantisse.size()-1; i++ )
    //    //std::cout << std::bitset<64>(pi.mantisse[ pi.mantisse.size() - i]);
    //std::cout << std::endl;
    pi.print();


}