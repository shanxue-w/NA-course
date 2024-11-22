// convert float to string
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <bitset>

float
Binary2Float(const std::string& binary)
{
    union 
    {
        float f;
        uint32_t i;
    } u;
    u.i = std::bitset<32>(binary).to_ulong();
    return u.f;
}

std::string
Float2Binary(const float f)
{
    union 
    {
        float f;
        uint32_t i;
    } u;
    u.f = f;
    return std::bitset<32>(u.i).to_string();
}

int main(void)
{
    std::cout << Float2Binary(1.0) << std::endl;
    std::cout << Float2Binary(0.9689124217106447) << std::endl;
    std::cout << Float2Binary(1.0-0.9689124217106447) << std::endl;
    return 0;
}

