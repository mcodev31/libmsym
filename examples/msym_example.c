//
//  msym_example.c
//  libmsym
//
//  Created by Marcus Johansson on 10/09/15.
//  Copyright (c) 2015 Marcus Johansson.
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

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
