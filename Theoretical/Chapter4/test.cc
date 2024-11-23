#include <iomanip>
#include <iostream>

int main(void) {
  float a = 1e-10;
  float b = 1e-5;
  float c = 1;

  float result = c;
  for (int i = 0; i < 100000; i++) {
    result += a;
  }

  std::cout << std::fixed << std::setprecision(20) << result << std::endl;

  result = 0.0;
  for (int i = 0; i < 100000; i++) {
    result += a;
  }
  result += c;
  std::cout << std::fixed << std::setprecision(20) << result << std::endl;

  // 十位小数
  std::cout.precision(10);
  std::cout << std::fixed << std::setprecision(20) << a + b + c << std::endl;
  std::cout.precision(10);
  std::cout << std::fixed << std::setprecision(20) << c + b + a << std::endl;
  return 0;
}