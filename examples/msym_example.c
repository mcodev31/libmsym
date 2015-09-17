
#include <stdio.h>
#include "example.h"

int main(int argc, const char * argv[]) {
    int ret = 1;
    if(argc == 2){
        ret = example(argv[1],NULL);
        fflush(stdout);
    } else {
        printf("usage msym_example <xyz-file>");
    }
    return ret;
}
