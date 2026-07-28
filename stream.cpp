#include <iostream>
#include <string>
#include <fstream>

int main(int argc, char **argv)
{
    using namespace std;
    std::ifstream fpfile("README.md");
    if (!fpfile)
        std::cerr << "Error opening file" << std::endl;
    else
        std::cout << "File opened successfully" << std::endl;
    fpfile.close();
}