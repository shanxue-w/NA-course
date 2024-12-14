#include "BInterpolate.hpp"
#include "PPInterpolate.hpp"
#include <fstream>     // For file I/O
#include <json/json.h> // For JSON handling

int
main(void)
{
    std::cout << "======================JSON======================" << std::endl;
    // Read the json file, control.json
    Json::Value  json;
    Json::Reader reader;

    // Open the file
    std::ifstream file("./src/control.json");
    if (!file.is_open())
    {
        std::cerr << "Failed to open file control.json" << std::endl;
        return 1;
    }

    bool parsingSuccessful = reader.parse(file, json);
    if (!parsingSuccessful)
    {
        std::cout << "Failed to parse configuration\n" << reader.getFormattedErrorMessages();
        return 1;
    }

    // Get the first element of the array (which is expected to be an object)
    const Json::Value &firstItem = json[0];

    // Now pass the first item (which is an object) to PPInterpolate
    PPInterpolate<3> pp(firstItem);
    BInterpolate<3>  bs(firstItem);
    double           maxerror_pp = 0;
    double           maxerror_bs = 0;
    for (double i = 0; i < 10; i += 0.1)
    {
        double error_pp = std::abs(pp(i) - 1.0);
        double error_bs = std::abs(bs(i) - 1.0);
        if (error_pp > maxerror_pp)
            maxerror_pp = error_pp;

        if (error_bs > maxerror_bs)
            maxerror_bs = error_bs;
    }

    std::cout << "Max error PP: " << maxerror_pp << std::endl;
    std::cout << "Max error BS: " << maxerror_bs << std::endl;
    std::cout << "======================JSON======================" << std::endl;
    return 0;
}
