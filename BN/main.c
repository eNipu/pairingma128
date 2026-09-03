#include "params.h"
#include "test.h"

int main(){
    init_parameters();
    set_parameters();
    print_parameters();

    test_ate_pairing();

    clear_parameters();
    return 0;
}
