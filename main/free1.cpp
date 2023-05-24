#include <functional>

template <class R, class... Args>
void hello<R(Args)>() {};

int main() {
    hello<double(double)>();
    return 0;
}